var geom = /* color: #98ff00 */ee.Geometry.Point([-76.21370770757636, 47.2501716861179]);
var collection = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR").filterBounds(geom);

function addL8EVI2(image) {
    // Bits 3 and 5 are cloud shadow and cloud, respectively.
    var cloudShadowBitMask = (1 << 3);
    var cloudsBitMask = (1 << 5);
    // Get the pixel QA band.
    var qa = image.select('pixel_qa');
    // Both flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
        .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
    var image = image.updateMask(mask);

    var evi2 = image.expression(
        '2.5 * ((0.0001 * NIR - 0.0001 * RED) / (0.0001 * NIR + 2.4 * 0.0001 * RED + 1))', {
          'NIR': image.select('B5'),
          'RED': image.select('B4'),
    });
    
    return image.addBands(evi2).set('system:time_start', image.get('system:time_start'));
}

var modcollection = ee.ImageCollection("MODIS/006/MCD43A4").filterBounds(geom).filterDate('2013-1-1', '2018-12-31');

function addMODISEVI2(image) {
    // Bitmask for BRDF_Albedo_Band_Mandatory_Quality_Band6
    var mask = image.select('BRDF_Albedo_Band_Mandatory_Quality_Band2').eq(0)
        .and(image.select('BRDF_Albedo_Band_Mandatory_Quality_Band2').eq(0))
    var image = image.updateMask(mask);
    var evi2 = image.expression(
        '2.5 * ((0.0001 * NIR - 0.0001 * RED) / (0.0001 * NIR + 2.4 * 0.0001 * RED + 1))', {
          'NIR': image.select('Nadir_Reflectance_Band2'),
          'RED': image.select('Nadir_Reflectance_Band1'),
    });    
    return image.addBands(evi2).set('system:time_start', image.get('system:time_start'));
}

// calculate EVI2 for the L8 data
var meanEVI2 = collection.map(addL8EVI2);


Map.addLayer(geom);
// print(ui.Chart.image.series(meanEVI2.select('constant'), geom, ee.Reducer.mean(), 30));