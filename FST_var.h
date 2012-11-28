
  ! FST variables 
  type fst_var
     character*4  shortname					! FST field name
     character*8  cdfname					! NetCDF field name
     character*64 longname					! description
     character*16 units						! unit of the field
  end type fst_var

  integer FST_VARS_COUNT
  parameter (FST_VARS_COUNT = 72+8+2)

  type (FST_var) :: vars(FST_VARS_COUNT)

data vars /  &
 FST_var('BR','BR' ,'Br concentration','kg/kg'), &
 FST_var('BRO','BRO','BrO concentration','kg/kg'), &
 FST_var('C001','GEM','GEM concentration','kg/kg'), & 
 FST_var('C024','RGM','RGM concentration','kg/kg'), & 
 FST_var('C003','TPM','TPM concentration','kg/kg'), & 
 FST_var('DWYT','DEP','Mercury total deposition','mcg/m2/year'), &
 FST_var('DEP','DEP','Mercury total deposition','mcg/m2/year'), &
 FST_var('DWTO','DEP','Mercury total deposition','mcg/m2'), &
 FST_var('DTOT','DEP','Mercury dry deposition','kg/m2'), &
 FST_var('WTOT','DEP','Mercury wet deposition','kg/m2'), &
 FST_var('AL','AL','Surface Albedo',''), & 
 FST_var('AP','AP','Planetary Albedo',''), & 
 FST_var('DD','DD','Isobaric Divergence','(/sec)'), &
 FST_var('GZ','Z','Geopotential Height','decametres'), &
 FST_var('ES','ES','Dew Point Depression','deg C'), &
 FST_var('TT','T','Air Temperature','deg C'), &
 FST_var('DS','DS','Divergence (sigma levels)','/sec'), &
 FST_var('DZ','DZ','Thickness','decametres'), &
 FST_var('FC','FC','Sensible Heat Flux','W/m2'), &
 FST_var('FQ','FQ','Momentum Flux','pascals'), &
 FST_var('FV','FV','Latent Heat Flux','W/m2'), &
 FST_var('GL','GL','Ice Cover',''), &
 FST_var( 'H', 'H','Boundary Layer Height','metres'), &
 FST_var('H1','H1','Boundary Layer Height','metres'), &
 FST_var('HR','HR','Relative Humidity',''), &
 FST_var('HS','HS','Potential Evaporation Fraction',''), &
 FST_var('HU','HU','Humidite specifique','kg/kg'), &
 FST_var('IZ','IZ','Infra-red Satellite Image',''), &
 FST_var('LA','LA','Geographical latitude',  'deg'  ), &
 FST_var('LO','LO','Geographical longitude', 'deg'  ), &
 FST_var('ME','ME','Topography','m'), &
 FST_var('MG','MG','Land-Sea Mask',''), &
 FST_var('MT','MT','Geopotential Topography','m2/s2'), &
 FST_var('MX','MX','Geopotential Topography','m2/s2'), &
 FST_var('MZ','MZ','Water Vapor Satellite Image',''), &
 FST_var('NB','NB','Low Cloud Cover',''), &
 FST_var('NE','NE','Snow Cover',''), &
 FST_var('NM','NM','Middle Cloud Cover',''), &
 FST_var('NH','NH','High Cloud Cover',''), &
 FST_var('NT','NT','Total Cloud Cover',''), &
 FST_var('P0','PS','Surface Pressure','millibars'), &
 FST_var('PN','PN','Sea Level Pressure','millibars'), &
 FST_var('PR','PR','Precipation Accumulation','metres'), &
 FST_var('PT','PT','Pressure at the Top of the Model','millibars'), &
 FST_var('PC','PC','Convective Precipation Accumulation','metres'), &
 FST_var('QG','QG','Ground humidity','kg/kg'), &
 FST_var('QQ','QQ','Absolute Vorticity','/sec'), &
 FST_var('QR','VO','Relative Vorticity','/sec'), &
 FST_var('QS','QS','Absolute Vorticity (sigma levels','/sec'), &
 FST_var('RN','RN','Stratiform Precipitation Accumulation','metres'), &
 FST_var('RR','RR','Liquid Precipitation Rate','m/s'), &
 FST_var('RC','RC','Convective Precipitation Rate','m/s'), &
 FST_var('RT','RT','Total Precipitation Rate','m/s'), &
 FST_var('SN','SN','Convective Precipitation Rate','metres'), &
 FST_var('SR','SR','Solid Precipitation Rate','m/s'), &
 FST_var('TI','TI','Infrared (IR) heating rate','K/day'), &
 FST_var('T2','T2','Visible (VIS) heating rate','K/day'), &
 FST_var('TG','TG','Ground Temperature','deg K'), &
 FST_var('TH','TH','Potential Temperature','deg K'), &
 FST_var('TM','TM','Sea Surface Temperature','deg K'), &
 FST_var('TP','TP','Deep Soil Temperature','deg C'), &
 FST_var('TR','TR','Deep Soil Temperature','deg C'), &
 FST_var('TS','TS','Surface Temperature','deg C'), &
 FST_var('TZ','TZ','Air Temperature','deg C'), &
 FST_var('UU','U','X Component of the Wind','m/s'), &
 FST_var('UV','UV','Wind Modulus','m/s'), &
 FST_var('VT','VT','Virtual Temperature','deg C'), &
 FST_var('VV','V','Y Component of the Wind','m/s'), &
 FST_var('VZ','VZ','Visibleble Satellite Image',''), &
 FST_var('WG','WG','Ground Humidity',''), &
 FST_var('WR','WR','Deep Soil Humidity',''), &
 FST_var('WW','OMEGA','Vertical Motion','pascals/s'), &
 FST_var('WS','WS','Vertical motion (sigma levels)','/sec'), &
 FST_var('Z0','Z0','Roughness Length','metres'), &
 FST_var('ZP','ZP','Log(Roughness Length)','log (metres)'), &
 FST_var('UT0','UT0','x component of velocity',''), &
 FST_var('VT0','VT0','y component of velocity',''), &
 FST_var('WT0','WT0','z component of velocity',''), &
 FST_var('TT0','TT0','temperature',''), &
 FST_var('TDT0','TDT0','total divergence (dpi* dot / dpi* + D )',''), &
 FST_var('FIT0','FIT0','geopotential',''), &
 FST_var('DV1','DV1','deposition velocity','m/s') /
