/**
 * This script borrows from and adds to code found in the Earth Engine Community tutorial
 * "Time Series Modeling" (https://developers.google.com/earth-engine/tutorials/community/time-series-modeling).

 * @license
 * Copyright 2022 The Google Earth Engine Community Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 *Significant changes include the addition of amplitude and phase calculations and 
 *the calculation of synthetic images based on the harmonic model (not currently part of the export).
 */

// An input geometry
var geometry = geometry3;
// Name for export file
var outname = 'Landsat_harmonic_regression';

// Import Landsat 8 image collection and define region of interest.
var landsat8Sr = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');

// The dependent variable we are modeling.
var dependent = 'NDVI'; // using NDVI here


// The number of cycles per year to model.
var harmonics = 3;

// Make a list of harmonic frequencies to model.  
// These also serve as band name suffixes.
var harmonicFrequencies = ee.List.sequence(1, harmonics);

// Function to get a sequence of band names for harmonic terms.
var getNames = function(base, list) {
  return ee.List(list).map(function(i) { 
    return ee.String(base).cat(ee.Number(i).int());
  });
};

// Construct lists of names for the harmonic terms (adding amplitude and phase)
var cosNames = getNames('cos_', harmonicFrequencies);
var sinNames = getNames('sin_', harmonicFrequencies);
var ampNames = getNames('amplitude_', harmonicFrequencies);
var phaseNames = getNames('phase_', harmonicFrequencies);

// Independent variables.
var independents = ee.List(['constant', 't']).cat(cosNames).cat(sinNames);

// Function to mask clouds in Landsat 8 imagery.
var maskL8sr = function(image) {
  // Bit 0 - Fill
  // Bit 1 - Dilated Cloud
  // Bit 2 - Cirrus
  // Bit 3 - Cloud
  // Bit 4 - Cloud Shadow
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // Apply the scaling factors to the appropriate bands.
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);

  // Replace the original bands with the scaled ones and apply the masks.
  return image.addBands(opticalBands, null, true)
    .addBands(thermalBands, null, true)
    .updateMask(qaMask)
    .updateMask(saturationMask);
};

// Function to add an NDVI band, the dependent variable.
var addNDVI = function(image) {
  return image
    .addBands(image.normalizedDifference(['SR_B5', 'SR_B4'])
    .rename('NDVI'))
    .float();
};

// Function to add a constant band.
var addConstant = function(image) {
  return image.addBands(ee.Image(1));
};

// Function to add a time band.
var addTime = function(image) {
  // Compute time in fractional years since the epoch.
  var date = image.date();
  var years = date.difference(ee.Date('1970-01-01'), 'year');
  var timeRadians = ee.Image(years.multiply(2 * Math.PI));
  return image.addBands(timeRadians.rename('t').float());
};

// Function to compute the specified number of harmonics
// and add them as bands. Assumes the time band is present.
var addHarmonics = function(freqs) {
  return function(image) {
    // Make an image of frequencies.
    var frequencies = ee.Image.constant(freqs);
    // This band should represent time in radians.
    var time = ee.Image(image).select('t');
    // Get the cosine terms.
    var cosines = time.multiply(frequencies).cos().rename(cosNames);
    // Get the sin terms.
    var sines = time.multiply(frequencies).sin().rename(sinNames);
    return image.addBands(cosines).addBands(sines);
  };
};




// Filter to the area of interest, mask clouds, add variables.
var harmonicLandsat = landsat8Sr
  .filterBounds(geometry)
  .filterDate('2013', '2024')
  .map(maskL8sr)
  .map(addNDVI)
  .map(addConstant)
  .map(addTime)
  .map(addHarmonics(harmonicFrequencies));
  
// The output of the regression reduction is a 4x1 array image.
var harmonicTrend = harmonicLandsat
  .select(independents.add(dependent)) //t, constant, cos_, sin_, NDVI
  .reduce(ee.Reducer.linearRegression(independents.length(), 1));

// Turn the array image into a multi-band image of coefficients.
var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([independents]);

// Compute fitted values.
var fittedHarmonic = harmonicLandsat.map(function(image) {
  return image.addBands(
    image.select(independents)
      .multiply(harmonicTrendCoefficients)
      .reduce('sum')
      .rename('fitted'));
});



// Add amplitutde
var harmonicTrendCoefficients = harmonicTrendCoefficients.addBands(harmonicTrendCoefficients.select('sin_.')
  .hypot(harmonicTrendCoefficients.select('cos_.'))
  // Add a scale factor for visualization.
  .multiply(5).rename(ampNames));

// Add phase 
var harmonicTrendCoefficients = harmonicTrendCoefficients.addBands(harmonicTrendCoefficients
  .select('sin_.').atan2(harmonicTrendCoefficients.select('cos_.'))
  .unitScale(-Math.PI, Math.PI).rename(phaseNames));

// Add Mean
var harmonicTrendCoefficients = harmonicTrendCoefficients
  .addBands(harmonicLandsat.select('NDVI').reduce('mean').rename('NDVI_mean'));

// https://www.sciencebase.gov/catalog/item/62169a62d34e4465d4091589
// Pull out the three bands for visualization.
var amplitude = harmonicTrendCoefficients.select('amplitude_1');
var phase = harmonicTrendCoefficients.select('phase_1');
var val = harmonicTrendCoefficients.select('NDVI_mean');
print(amplitude, "amplitude");
// Turn the HSV data into an RGB image and add it to the map.
var seasonality = ee.Image.cat(phase, amplitude, val).hsvToRgb();


// Create composites
var stepList = ee.List.sequence(2013,2023);
var filterCollection = stepList.map(function(year){
  var startDate = ee.Date.fromYMD(year,6,1);
  var endDate = ee.Date.fromYMD(year,9,30);
  var composite_i = fittedHarmonic.filterDate(startDate, endDate)
                        .median()
                        .set('system:time_start', startDate);
    
  // add a 'valid' property of 1 if image has bands and 0 if image has no bands
  // it works because a populated list is evaluated as true and an empty list is evaluated as false
  composite_i = ee.Algorithms.If(composite_i.bandNames(), composite_i.set('valid', 1), composite_i.set('valid', 0));
    
    return composite_i;
});

// make an image collection from list of images
var yearlyComposites = ee.ImageCollection(filterCollection);

// filter only valid images
var yearlyComposites = yearlyComposites.filterMetadata('valid', 'equals', 1);
print(yearlyComposites, "yearlyComposites");
var yearlyNDVI = yearlyComposites.select('fitted');

var outImage = yearlyNDVI.toBands()
        .select(["0_fitted", "1_fitted", "2_fitted", "3_fitted", "4_fitted", "5_fitted", "6_fitted", "7_fitted", "8_fitted", "9_fitted", "10_fitted"],
        ["2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"])
        .multiply(10000)
        .toShort(); // cast to 16 bit signed integer

// Prepare NDVI fitted output (synthetic) images
// Currently turned off..
// var outImage = yearlyNDVI.toBands()
//         .multiply(10000)
//         .toShort(); // cast to 16 bit signed integer
// print(outImage, 'NDVI Collection to image');

Map.addLayer(outImage.select('2021'), {min: -10000, max: 10000, palette: ['red', 'brown', 'yellow', 'green']}, 'yearly NDVI 2021');   
Map.addLayer(outImage.select('2013'), {min: -10000, max: 10000, palette: ['red', 'brown', 'yellow', 'green']}, 'yearly NDVI 2013');
// Display the ROI and NDVI composite on the map.
Map.addLayer(seasonality, {}, 'Seasonality Harmonic 1');

Map.addLayer(harmonicLandsat,
  {bands: 'NDVI', min: 0.1, max: 0.9, palette: ['white', 'green']},
  'NDVI Mosaic');
// Display amplitude
Map.addLayer(harmonicTrendCoefficients,
  {bands: 'amplitude_1', min: 0.1, max: 1.75, palette: ['white', 'green']},
  'Harmonic 1 Amplitude');
// Map.addLayer(harmonicTrendCoefficients,
//   {bands: 'amplitude_6', min: 0.1, max: 1.75, palette: ['white', 'green']},
//   'Harmonic 6 Amplitude');
  
  
Map.addLayer(fittedHarmonic.first(),
  {bands: 'fitted', min: 0.1, max: 1.75, palette: ['white', 'green']},
  'fitted harmonic');
print(fittedHarmonic, "fittedHarmonic");
print(harmonicTrend, "harmonicTrend");
print(harmonicTrendCoefficients, "harmonicTrendCoeff");
print(harmonicLandsat, "harmonicLandsat");
// Map.addLayer(fittedHarmonic, {bands: 'fitted'});

var output = harmonicTrendCoefficients.multiply(10000).toShort();
var output = output.clip(geometry);
// Export harmonic coefficients sin, cos, amp, phase. Not exporting synthetic images.
Export.image.toDrive({
  folder: 'GEE_Landsat_Exports',
  image: output,
  description: outname,
  fileNamePrefix: outname, 
  region: geometry,
  crs: "EPSG:5070",
  scale: 30,
  fileFormat: 'GeoTIFF',
  // formatOptions: {
  //   noData: noDataVal,
  // },
  // set the x and y translation here based on our template raster for current project.
  crsTransform: [30.0, 0.0, 4.004536, 0.0, -30.0, -13.497535],
  maxPixels: 1e13
});


