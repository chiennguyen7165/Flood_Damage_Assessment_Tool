/*===========================================================================================

  *******************************************************************************************
                                  LƯA CHỌN KHU VỰC NGHIÊN CỨU   

   Mặc định là tỉnh Quảng Trị, shape file đã được import

  *******************************************************************************************
                                       LỰA CHỌN THỜI GIAN ĐÁNH GIÁ
  /*
   Lựa chọn khoảng thời gian trước và sau khi lũ lụt xảy ra. Nên để khoảng thời gian đủ dài để thu
   thập ảnh Sentinel-1. Khuyến nghị từ 15 ngày trở nên
   */

// Thời điểm bắt đầu và kết thúc của giai đoạn trước khi xảy ra lũ lụt
var before_start= '2021-06-06';
var before_end='2021-10-06';

// Thời điểm bắt đầu và kết thúc của giai đoạn sau khi xảy ra lũ lụt
var after_start='2021-10-06';
var after_end='2022-02-06';

/********************************************************************************************
                           THIẾT LẬP CÁC THÔNG SỐ CHẠY (Nên giữ nguyên)*/

// Lựa chọn chiều phân cực của ảnh VV; HH; hoặc VH
var polarization = "VH"; 

// Chế độ chiếu của ảnh SAR. ASCENDING và DESCENDING.
var pass_direction = "DESCENDING"; 

// Ngưỡng phân biệt giữa pixel nước và các pixel khác. Có thể hiệu chỉnh hoặc để nguyên
var difference_threshold = 1.25; 

//*********************************************************************************************

//---------------------------------- Lấy giá trị đầu vào của người dùng ------------------------------//
//--------------------------------------- QUÁ TRÌNH TIỀN XỬ LÝ ẢNH --------------------------//

// Thay đổi lại tên biến chứa khu vực nghiên cứu
var aoi = ee.FeatureCollection(geometry);

// Load và lọc ảnh Sentinel-1 theo các biến setup
var collection= ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.eq('instrumentMode','IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', polarization))
  .filter(ee.Filter.eq('orbitProperties_pass',pass_direction)) 
  .filter(ee.Filter.eq('resolution_meters',10))
  .filterBounds(aoi)
  .select(polarization);
  
// Lọc ảnh theo dữ liệu thời gian
// biến chứa các ảnh trước lũ
var before_collection = collection.filterDate(before_start, before_end);

// biến chứa các ảnh sau lũ
var after_collection = collection.filterDate(after_start,after_end);

// In kết quả lên màn hình

      // Lọc kết quả ngày có ảnh
      function dates(imgcol){
        var range = imgcol.reduceColumns(ee.Reducer.minMax(), ["system:time_start"]);
        var printed = ee.String('Từ ')
          .cat(ee.Date(range.get('min')).format('dd-MM-YYYY'))
          .cat(' đến ')
          .cat(ee.Date(range.get('max')).format('dd-MM-YYYY'));
        return printed;
      }
      // In kết quả thông tin ảnh lọc được trước lũ lên màn hình console
      var before_count = before_collection.size();
      print(ee.String('Số ảnh lọc được: Trước lũ ').cat('(').cat(before_count).cat(')'),
        dates(before_collection), before_collection);
      
      // In kết quả thông tin ảnh lọc được sau lũ lên màn hình console
      var after_count = before_collection.size();
      print(ee.String('Số ảnh lọc được: Sau lũ ').cat('(').cat(after_count).cat(')'),
        dates(after_collection), after_collection);

// Tiến hành gộp các bộ sưu tập ảnh trước và sau lũ lại thành 1 ảnh tương ứng với mỗi giai đoạn
var before = before_collection.mosaic().clip(aoi);
var after = after_collection.mosaic().clip(aoi);

// Lọc các pixel gây nhiễu ảnh
var smoothing_radius = 50;
var before_filtered = before.focal_mean(smoothing_radius, 'circle', 'meters');
var after_filtered = after.focal_mean(smoothing_radius, 'circle', 'meters');


//------------------------------- TÍNH TOÁN DIỆN TÍCH LŨ -------------------------------//

// Tính toán sự khác nhau giữa hai ảnh trước và sau lũ
var difference = after_filtered.divide(before_filtered);

// Lọc các pixel là nước bằng cách sử dụng ngưỡng phân cực
var threshold = difference_threshold;
var difference_binary = difference.gt(threshold);

// Hiệu chỉnh tiếp kết quả bằng các bộ sưu tập ảnh khác
      
      // Lấy BST JRC layer về nước mặt vệ tinh để lọc các pixel là nước dài hạn
      // (các pixel đó là nước được hơn 10 tháng trong 1 năm)
      var swater = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').select('seasonality');
      var swater_mask = swater.gte(10).updateMask(swater.gte(10));
      
      //Lọc kết quả bằng điều kiện trên (water > 10 mo/yr) các pixel 
      // thỏa mãn điều kiện là các pixel có giá trị  = 0
      var flooded_mask = difference_binary.where(swater_mask,0);
      // Rút ra kết quả cuối cùng thể hiện diện tích ngập lụt
      var flooded = flooded_mask.updateMask(flooded_mask);
      
      // Lọc nhiễu cho kết quả
      var connections = flooded.connectedPixelCount();    
      var flooded = flooded.updateMask(connections.gte(8));
      
      // Che dấu các khu vực có độ dốc hơn 5 phần trăm bằng bộ sưu tập về độ cao kỹ thuật DEM
      var DEM = ee.Image('WWF/HydroSHEDS/03VFDEM');
      var terrain = ee.Algorithms.Terrain(DEM);
      var slope = terrain.select('slope');
      var flooded = flooded.updateMask(slope.lt(5));

// Tính toán diện tích ngập lụt
// Tạo lớp raster chứa thông tin về diện tích ngập lục của từng pixel
var flood_pixelarea = flooded.select(polarization)
  .multiply(ee.Image.pixelArea());

// Tính tổng diện tích của các pixel bị ngập
var flood_stats = flood_pixelarea.reduceRegion({
  reducer: ee.Reducer.sum(),              
  geometry: aoi,
  scale: 10,
  bestEffort: true
  });

// Chuyển đổi diện tích ngập lụt thành 
var flood_area_ha = flood_stats
  .getNumber(polarization)
  .divide(10000)
  .round(); 


//------------------------------  ĐÁNH GIÁ THIỆT HẠI  ----------------------------------//

//----------------------------- Số người bị ảnh hưởng ----------------------------//

// Load bộ sưu tập JRC Global Human Settlement Popluation Density layer: mật độ dân số toàn cầu
// Độ phân giải: 250. Số người trên mỗi ô được cung cấp.
var population_count = ee.Image('JRC/GHSL/P2016/POP_GPW_GLOBE_V1/2015').clip(aoi);

// Tạo lớp chứa thông tin về số người dân
var GHSLprojection = population_count.projection();

// Tính toán lại diện tích ngập lụt lọc theo thông tin của biến người dân vừa tạo
var flooded_res1 = flooded
    .reproject({
    crs: GHSLprojection
  });

// Tạo một lớp raster mới thể hiện dân số ảnh hưởng mà chỉ sử dụng lớp lũ mới được tính toán lại
var population_exposed = population_count
  .updateMask(flooded_res1)
  .updateMask(population_count);

// Tính tổng giá trị số người bị ảnh hưởng của các pixcel 
var stats = population_exposed.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: aoi,
  scale: 250,
  maxPixels:1e9 
});

// Chuyển đổi kết quả tính toán sang dạng số nguyên
var number_pp_exposed = stats.getNumber('population_count').round();

//----------------------------- DIỆN TÍCH ĐẤT NÔNG NGHIỆP BỊ ẢNH HƯỞNG ----------------------------//

// Load bộ sưu tập ảnh MODIS về lớp phủ đất mặt toàn cầu
var LC = ee.ImageCollection('MODIS/006/MCD12Q1')
  .filterDate('2014-01-01',after_end)
  .sort('system:index',false)
  .select("LC_Type1")
  .first()
  .clip(aoi);

// Tạo biến chứa thông tin chiết tách toàn bộ các pixel là cây trồng
var cropmask = LC
  .eq(12)
  .or(LC.eq(14))
var cropland = LC
  .updateMask(cropmask)
  
// Chiết tách ảnh MODIS bằng biến điều kiện vừa tạo
var MODISprojection = LC.projection();

// Tính toán lại diện tích ngập lụt theo kết quả chiết tách ảnh MODIS
var flooded_res = flooded
    .reproject({
    crs: MODISprojection
  });

// Tính toán diện tích đất trồng trọt bị ảnh hưởng bằng cách sử dụng lớp lũ được tính toán lại
var cropland_affected = flooded_res
  .updateMask(cropland)

// Lấy diện tích các pixel là cây trồng bị ảnh hưởng
var crop_pixelarea = cropland_affected
  .multiply(ee.Image.pixelArea()); //tính toán diện tích của từng pixel

// Tính tổng giá trị của các pixel
var crop_stats = crop_pixelarea.reduceRegion({
  reducer: ee.Reducer.sum(),              
  geometry: aoi,
  scale: 500,
  maxPixels: 1e9
  });
  
// Chuyển đổi đơn vị sang hecta
var crop_area_ha = crop_stats
  .getNumber(polarization)
  .divide(10000)
  .round();


//------------------------------  DISPLAY PRODUCTS  ----------------------------------//

// Before and after flood SAR mosaic
Map.centerObject(aoi,8);
Map.addLayer(before_filtered, {min:-25,max:0}, 'Trước lũ',0);
Map.addLayer(after_filtered, {min:-25,max:0}, 'Sau lũ',1);

// Difference layer
Map.addLayer(difference,{min:0,max:2},"Lớp nước khác nhau",0);

// Flooded areas
Map.addLayer(flooded,{palette:"0000FF"},'Diện tích ngập');

// Population Density
var populationCountVis = {
  min: 0,
  max: 200.0,
  palette: ['060606','337663','337663','ffffff'],
};
Map.addLayer(population_count, populationCountVis, 'Mật độ dân số',0);

// Exposed Population
var populationExposedVis = {
  min: 0,
  max: 200.0,
  palette: ['yellow', 'orange', 'red'],
};
Map.addLayer(population_exposed, populationExposedVis, 'Mật độ dân số bị ảnh hưởng');

// MODIS Land Cover
var LCVis = {
  min: 1.0,
  max: 17.0,
  palette: [
    '05450a', '086a10', '54a708', '78d203', '009900', 'c6b044', 'dcd159',
    'dade48', 'fbff13', 'b6ff05', '27ff87', 'c24f44', 'a5a5a5', 'ff6d4c',
    '69fff8', 'f9ffa4', '1c0dff'
  ],
};
Map.addLayer(LC, LCVis, 'Lớp phủ đất',0);

// Cropland
var croplandVis = {
  min: 0,
  max: 14.0,
  palette: ['30b21c'],
};
Map.addLayer(cropland, croplandVis, 'Đất nông nghiệp',0)

// Affected cropland
Map.addLayer(cropland_affected, croplandVis, 'Đất nông nghiệp bị ảnh hưởng'); 



//------------------------------------- EXPORTS ------------------------------------//
// Export flooded area as TIFF file 
Export.image.toDrive({
  image: flooded, 
  description: 'Flood_extent_raster',
  fileNamePrefix: 'flooded',
  region: aoi, 
  maxPixels: 1e10
});

// Export flooded area as shapefile (for further analysis in e.g. QGIS)
// Convert flood raster to polygons
var flooded_vec = flooded.reduceToVectors({
  scale: 10,
  geometryType:'polygon',
  geometry: aoi,
  eightConnected: false,
  bestEffort:true,
  tileScale:2,
});

// Export flood polygons as shape-file
Export.table.toDrive({
  collection:flooded_vec,
  description:'Flood_extent_vector',
  fileFormat:'SHP',
  fileNamePrefix:'flooded_vec'
});

// Export auxcillary data as shp?
// Exposed population density
Export.image.toDrive({
  image:population_exposed,
  description:'Exposed_Populuation',
  scale: 250,
  fileNamePrefix:'population_exposed',
  region: aoi,
  maxPixels:1e10
});

//---------------------------------- MAP PRODUCTION --------------------------------//

//-------------------------- Display the results on the map -----------------------//

// set position of panel where the results will be displayed 
var results = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px',
    width: '350px'
  }
});

//Prepare the visualtization parameters of the labels 
var textVis = {
  'margin':'0px 8px 2px 0px',
  'fontWeight':'bold'
  };
var numberVIS = {
  'margin':'0px 0px 15px 0px', 
  'color':'bf0f19',
  'fontWeight':'bold'
  };
var subTextVis = {
  'margin':'0px 0px 2px 0px',
  'fontSize':'12px',
  'color':'grey'
  };

var titleTextVis = {
  'margin':'0px 0px 15px 0px',
  'fontSize': '18px', 
  'font-weight':'', 
  'color': '3333ff'
  };

// Create lables of the results 
// Titel and time period
var title = ui.Label('Kết quả', titleTextVis);
var text1 = ui.Label('Tình hình lũ lụt giai đoạn:',textVis);
var number1 = ui.Label(after_start.concat(" - ",after_end),numberVIS);


// Estimated flood extent 
var text2 = ui.Label('Ước tính diện tích ngập lụt:',textVis);
var text2_2 = ui.Label('Đang ước tính...',subTextVis);
dates(after_collection).evaluate(function(val){text2_2.setValue('theo dữ liệu ảnh Senintel-1 '+val)});
var number2 = ui.Label('Đang ước tính...',numberVIS); 
flood_area_ha.evaluate(function(val){number2.setValue(val+' hecta')}),numberVIS;

// Estimated number of exposed people
var text3 = ui.Label('Ước tính số người bị ảnh hưởng bởi lũ: ',textVis);
var text3_2 = ui.Label('theo dữ liệu GHSL 2015 (250m)',subTextVis);
var number3 = ui.Label('Đang ước tính...',numberVIS);
number_pp_exposed.evaluate(function(val){number3.setValue(val)}),numberVIS;

// Estimated area of affected cropland 
var MODIS_date = ee.String(LC.get('system:index')).slice(0,4);
var text4 = ui.Label('Ước tính diện tích đất nông nghiệp bị ảnh hưởng:',textVis);
var text4_2 = ui.Label('Đang ước tính', subTextVis)
MODIS_date.evaluate(function(val){text4_2.setValue('theo dữ liệu ảnh MODIS '+val +' (500m)')}), subTextVis;
var number4 = ui.Label('Đang ước tính...',numberVIS);
crop_area_ha.evaluate(function(val){number4.setValue(val+' hecta')}),numberVIS;

// Giới thiệu
var thesis = 'Sản phẩm của đề tài: Ứng dụng công nghệ điện toán đám mây Google Earth Engine đánh giá thiệt hại sau đợt lũ năm 2020 tại tỉnh Quảng Trị'
var text6 = ui.Label(thesis,subTextVis)

// Tác gia
var text7 = ui.Label('Học viên: Nguyễn Quang Chiến', subTextVis)

// Add the labels to the panel 
results.add(ui.Panel([
        title,
        text1,
        number1,
        text2,
        text2_2,
        number2,
        text3,
        text3_2,
        number3,
        text4,
        text4_2,
        number4,
        text6,
        text7]
      ));

// Add the panel to the map 
Map.add(results);

//----------------------------- Display legend on the map --------------------------//

// Tạo chú thích
var legend = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px 15px',
  }
});
 
// Create legend title
var legendTitle = ui.Label('Chú thích',titleTextVis);
 
// Add the title to the panel
legend.add(legendTitle);
 
// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
 
      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
 
      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
 
//  Palette with the colors
var palette =['#0000FF', '#30b21c'];
 
// name of the legend
var names = ['Diện tích ngập lụt','Diện tích đất nông nghiệp'];
 
// Add color and and names
for (var i = 0; i < 2; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  

// Create second legend title to display exposed population density
var legendTitle2 = ui.Label({
value: 'Mật độ dân số',
style: {
fontWeight: 'bold',
fontSize: '15px',
margin: '10px 0 0 0',
padding: '0'
}
});

// Add second title to the panel
legend.add(legendTitle2);

// create the legend image
var lon = ee.Image.pixelLonLat().select('latitude');
var gradient = lon.multiply((populationExposedVis.max-populationExposedVis.min)/100.0).add(populationExposedVis.min);
var legendImage = gradient.visualize(populationExposedVis);
 
// create text on top of legend
var panel = ui.Panel({
widgets: [
ui.Label('> '.concat(populationExposedVis['max']))
],
});
 
legend.add(panel);
 
// create thumbnail from the image
var thumbnail = ui.Thumbnail({
image: legendImage,
params: {bbox:'0,0,10,100', dimensions:'10x50'},
style: {padding: '1px', position: 'bottom-center'}
});
 
// add the thumbnail to the legend
legend.add(thumbnail);
 
// create text on top of legend
var panel = ui.Panel({
widgets: [
ui.Label(populationExposedVis['min'])
],
});
 
legend.add(panel);
 
// add legend to map (alternatively you can also print the legend to the console)
Map.add(legend);
