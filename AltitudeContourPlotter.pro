PRO DataCubeContourPlotter
;;Save pictures location
picsavloc = "C:\Users\rshar\Pictures"
;;Save output location
outputsavloc = 'C:\Users\rshar\IDLWorkspace83\Default\Current Data\April 2017\2D Plot Test Area\Processed Data\'
  
;;Altitude Version


file = dialog_pickfile(path = 'C:\Users\rshar\IDLWorkspace83\Default\Current Data\April 2017\2D Plot Test Area')
restore, file

Variables = ['Temperature','Co2', 'Co','N2','O','Co-Co2 Ratio', 'N2-Co2 Ratio','O-Co2 Ratio']
VarNameLength = n_elements(Variables)

print, '1: Altitude'
print, '2: Longitude'
print, '3: Latitude'

ask1 = float(ask('the number corresponding to the variable that you would like to be the Y coordinate '))
ask2 = float(ask('the number corresponding to the variable that you would like to be the X coordinate '))
ask3 = float(ask('the number corresponding to the variable that you would like to be the Z coordinate '))

for v = 0,(VarNameLength-1) do begin
  print, strtrim(v+1,1), ': ', Variables[v]
endfor

ask4 = float(ask('the corresponding number for the Variable you would like to plot'))

miny = float(ask('the minimum altitude you would like to look at'))
maxy = float(ask('the maximum altitude you would like to look at'))

if ask3 eq 3 then begin
  minz = float(ask('the minimum Latitude in your range'))
  maxz = float(ask('the maximum Latitude in your range'))
  namez = 'Lat'
  namex = 'Lon'
endif

if ask3 eq 2 then begin
  minz = float(ask('the minimum Longitude in your range'))
  maxz = float(ask('the maximum Longitude in your range'))
  namez = 'Lon'
  namex = 'Lat'
endif

name = string(ask('a name for this run if you would like to save the processed data. Otherwise, enter "n" '))

;Correcting user input to match grid for altitude
roundingfactor = 1.25
divminy = miny/roundingfactor
trueminy = ceil(divminy)*roundingfactor

if abs(trueminy mod 2) eq 0 then begin
  trueminy = trueminy - roundingfactor
endif

divmaxy = maxy/roundingfactor
truemaxy = ceil(divmaxy)*roundingfactor

if abs(truemaxy mod 2) eq 0 then begin
  truemaxy = truemaxy - roundingfactor
endif

print, 'Altitude range: ',  trueminy,' to ', truemaxy


miny = 1*trueminy
maxy = 1*truemaxy

;Correcting user input to match grid for latitude
if ask3 eq 3 then begin
  roundingfactor2 = 2.5
  divminz = minz/roundingfactor2
  trueminz = ceil(divminz)*roundingfactor2

  if abs(trueminz mod 2) eq 1 then begin
    trueminz = trueminz - roundingfactor2
  endif

  divmaxz = maxz/roundingfactor2
  truemaxz = ceil(divmaxz)*roundingfactor2

  if abs(truemaxz mod 2) eq 1 then begin
    truemaxz = truemaxz - roundingfactor2
  endif

  print, 'Latitude range: ',  trueminz,' to ', truemaxz

  minz = 1*trueminz
  maxz = 1*truemaxz  
endif

;Correcting user input to match grid for longitude
if ask3 eq 2 then begin
  roundingfactor2 = 2.5
  divminz = minz/roundingfactor2
  trueminz = ceil(divminz)*roundingfactor2

  if abs(trueminz mod 2) eq 0 then begin
    trueminz = trueminz + roundingfactor2
  endif

  divmaxz = maxz/roundingfactor2
  truemaxz = ceil(divmaxz)*roundingfactor2

  if abs(truemaxz mod 2) eq 0 then begin
    truemaxz = truemaxz - roundingfactor2
  endif

  print, 'Longitude range: ',  trueminz,' to ', truemaxz

  minz = 1*trueminz
  maxz = 1*truemaxz
endif


if ask2 eq 2 then begin
  minx = 2.5
  maxx = 357.50
endif

if ask2 eq 3 then begin
  minx = -87.50
  maxx = 87.50
endif


minalt = 1.25
binalt = 2.5
maxalt = 298.75

minlon = 2.5
binlon = 5.0
maxlon = 357.50

minlat = -87.50
binlat = 5.0
maxlat = 87.50

if ask3 eq 3 then begin
  xcoordsteps = (maxx - minx)/binlon + 1
endif

if ask3 eq 2 then begin
  xcoordsteps = (maxx - minx)/binlat + 1
endif


if ask3 eq 3 then begin
  lowestlatitude = minz
  Zdim = lowestlatitude*1
  Zdimp = 1*Zdim
  highestlatitude = maxz
endif

if ask3 eq 2 then begin
  lowestlongitude = minz
  Zdim = lowestlongitude*1
  Zdimp = 1*Zdim
  highestlongitude = maxz
endif

ycoordsteps = (maxy - miny)/binalt + 1




X = []
factor = 0
top = MAKE_ARRAY(xcoordsteps, ycoordsteps, /float, VALUE = 0.0)
top2 = MAKE_ARRAY(xcoordsteps, ycoordsteps, /float, VALUE = 0.0)

if ask3 eq 3 then begin
  
  latitude_slices = abs((lowestlatitude - highestlatitude)/binlat)
  for n = 0,(latitude_slices) do begin
    
    ;print, 'Averaging slice at', Zdim
  
    ;;Setup the Y Coordinate according to user input
    if ask1 eq 1 then begin
      yrange1 = [minalt,maxalt]
      ybinsize = binalt
      altminloc = fix((miny - minalt)/binalt)
      altmaxloc = fix((maxy - minalt)/binalt)
    endif
  
    if ask1 eq 2 then begin
      yrange1 = [minlon,maxlon]
      ybinsize = binlon
      lonminloc = fix((miny - minlon)/binlon)
      lonmaxloc = fix((maxy - minlon)/binlon)
    endif
  
    if ask1 eq 3 then begin
      yrange1 = [minlat,maxlat]
      ybinsize = binlat
      latminloc = fix((miny - minlat)/binlat)
      latmaxloc = fix((maxy - minlat)/binlat)
    endif
  
  
    ;;Setup the X Coordinate according to user input
    if ask2 eq 1 then begin
      xrange1 = [minalt,maxalt]
      xbinsize = binalt
      altminloc = fix((minx - minalt)/binalt)
      altmaxloc = fix((maxx - minalt)/binalt)
    endif
  
    if ask2 eq 2 then begin
      xrange1 = [minlon,maxlon]
      xbinsize = binlon
      lonminloc = fix((minx - minlon)/binlon)
      lonmaxloc = fix((maxx - minlon)/binlon)
    endif
  
    if ask2 eq 3 then begin
      xrange1 = [minlat,maxlat]
      xbinsize = binlat
      latminloc = fix((minx - minlat)/binlat)
      latmaxloc = fix((maxx - minlat)/binlat)
    endif
  
    ;;Setup the Z coordinate according to user input
    if ask3 eq 1 then begin
      altminloc = fix(Zdim - minalt)/binalt
      altmaxloc = fix(Zdim - minalt)/binalt
    endif
  
    if ask3 eq 2 then begin
      lonminloc = fix(Zdim - minlon)/binlon
      lonmaxloc = fix(Zdim - minlon)/binlon
    endif
  
    if ask3 eq 3 then begin
      latminloc = fix(Zdim - minlat)/binlat
      latmaxloc = fix(Zdim - minlat)/binlat
    endif
  
  
    ; Set up variables for the plot. For now, using Altitude as the Y coordinate and Longitude as the X coordinate
    if ask4 eq 1 then begin
      densityslice_i = reform(temperature.Tn[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      varname = 'Temperature'
    endif
    
    if ask4 eq 2 then begin
      varname = 'Co2'
      densityslice_i = reform(ndensitys.co2[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    endif
    
    if ask4 eq 3 then begin
      varname = 'Co'
      densityslice_i = reform(ndensitys.Co[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    endif
    
    if ask4 eq 4 then begin
      varname = 'N2'
      densityslice_i = reform(ndensitys.N2[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    endif
    
    if ask4 eq 5 then begin
      varname = 'Oxygen'
      densityslice_i = reform(ndensitys.O[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    endif
    
    if ((ask4 ge 6) and (ask4 le 8)) then begin
      varname = Variables[(ask4-1)]
      densityslice_i = reform(ndensitys.co2[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      if ask4 eq 6 then begin
        densityslice2_i = reform(ndensitys.co[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      endif
      if ask4 eq 7 then begin
        densityslice2_i = reform(ndensitys.N2[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      endif
      if ask4 eq 8 then begin
        densityslice2_i = reform(ndensitys.O[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      endif
    endif
    factor_i = cos(Zdim*!dtor)
    
    top_i = factor_i*densityslice_i
    
    if (ask4 ge 6) and (ask4 le 8) then begin
      top2_i = factor_i*densityslice2_i
    endif
    
    factor = factor + factor_i
    
    top = top + top_i
    if ((ask4 ge 6) and (ask4 le 8)) then begin
      top2 = top2 + top2_i
    endif
    
    Zdim = Zdim + binlat
  endfor
endif

;;Longitude
if ask3 eq 2 then begin

  longitude_slices = abs((lowestlongitude - highestlongitude)/binlon)
  for n = 0,(longitude_slices) do begin

    ;print, 'Averaging slice at', Zdim

    ;;Setup the Y Coordinate according to user input
    if ask1 eq 1 then begin
      yrange1 = [minalt,maxalt]
      ybinsize = binalt
      altminloc = fix((miny - minalt)/binalt)
      altmaxloc = fix((maxy - minalt)/binalt)
    endif

    if ask1 eq 2 then begin
      yrange1 = [minlon,maxlon]
      ybinsize = binlon
      lonminloc = fix((miny - minlon)/binlon)
      lonmaxloc = fix((maxy - minlon)/binlon)
    endif

    if ask1 eq 3 then begin
      yrange1 = [minlat,maxlat]
      ybinsize = binlat
      latminloc = fix((miny - minlat)/binlat)
      latmaxloc = fix((maxy - minlat)/binlat)
    endif


    ;;Setup the X Coordinate according to user input
    if ask2 eq 1 then begin
      xrange1 = [minalt,maxalt]
      xbinsize = binalt
      altminloc = fix((minx - minalt)/binalt)
      altmaxloc = fix((maxx - minalt)/binalt)
    endif

    if ask2 eq 2 then begin
      xrange1 = [minlon,maxlon]
      xbinsize = binlon
      lonminloc = fix((minx - minlon)/binlon)
      lonmaxloc = fix((maxx - minlon)/binlon)
    endif

    if ask2 eq 3 then begin
      xrange1 = [minlat,maxlat]
      xbinsize = binlat
      latminloc = fix((minx - minlat)/binlat)
      latmaxloc = fix((maxx - minlat)/binlat)
    endif

    ;;Setup the Z coordinate according to user input
    if ask3 eq 1 then begin
      altminloc = fix(Zdim - minalt)/binalt
      altmaxloc = fix(Zdim - minalt)/binalt
    endif

    if ask3 eq 2 then begin
      lonminloc = fix(Zdim - minlon)/binlon
      lonmaxloc = fix(Zdim - minlon)/binlon
    endif

    if ask3 eq 3 then begin
      latminloc = fix(Zdim - minlat)/binlat
      latmaxloc = fix(Zdim - minlat)/binlat
    endif


    ; Set up variables for the plot. For now, using Altitude as the Y coordinate and Longitude as the X coordinate
    if ask4 eq 1 then begin
      densityslice_i = reform(temperature.Tn[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      varname = 'Temperature'
    endif

    if ask4 eq 2 then begin
      varname = 'Co2'
      densityslice_i = reform(ndensitys.co2[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    endif

    if ask4 eq 3 then begin
      varname = 'Co'
      densityslice_i = reform(ndensitys.Co[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    endif

    if ask4 eq 4 then begin
      varname = 'N2'
      densityslice_i = reform(ndensitys.N2[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    endif

    if ask4 eq 5 then begin
      varname = 'Oxygen'
      densityslice_i = reform(ndensitys.O[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    endif

    if ((ask4 ge 6) and (ask4 le 8)) then begin
      varname = Variables[(ask4-1)]
      densityslice_i = reform(ndensitys.co2[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      if ask4 eq 6 then begin
        densityslice2_i = reform(ndensitys.co[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      endif
      if ask4 eq 7 then begin
        densityslice2_i = reform(ndensitys.N2[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      endif
      if ask4 eq 8 then begin
        densityslice2_i = reform(ndensitys.O[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
      endif
    endif
    factor_i = 1

    top_i = factor_i*densityslice_i

    if (ask4 ge 6) and (ask4 le 8) then begin
      top2_i = factor_i*densityslice2_i
    endif

    factor = factor + factor_i

    top = top + top_i
    if ((ask4 ge 6) and (ask4 le 8)) then begin
      top2 = top2 + top2_i
    endif

    Zdim = Zdim + binlat
  endfor
endif

   dens1 = top/factor
   dens2 = top2/factor
  
  if ((ask4 ge 6) and (ask4 le 8)) then begin
   X = dens2/dens1
   print, 'Ratio'
  endif else begin
   X = dens1
   print, 'No Ratio'
  endelse
   
   
   Xmax = max(X)
   Xmin = min(X)
   
   alts = meta.altitude[altminloc:altmaxloc]
   lons = meta.longitude[lonminloc:lonmaxloc]
   lats = meta.latitude[latminloc:latmaxloc]
   
   minlonaxis = min(lons)
   maxlonaxis = max(lons)
   
   minaltaxis = min(alts)
   maxaltaxis = max(alts)
   
   minlataxis = min(lats)
   maxlataxis = max(lats)
   
   if ask3 eq 3 then begin
     xrange1 = [minlonaxis,maxlonaxis]
     yrange1 = [minaltaxis,maxaltaxis]
   endif
   
   if ask3 eq 2 then begin
     xrange1 = [minlataxis,maxlataxis]
     yrange1 = [minaltaxis,maxaltaxis]
   endif
   
   ; Create the density plot by binning the data into a 2D histogram.
   if ((ask4 eq 1)) then begin
    density = X
   endif else begin
    density = log(X)
   endelse
    
;   if ((ask4 ge 6) and (ask4 le 8)) then begin
;    RatioEqual1 = MAKE_ARRAY(xcoordsteps, ycoordsteps, /float, VALUE = log(-1))
;    RatioEqual1[WHERE(((density gt -.05) and (density lt .05)), /NULL)] = 1 
;   endif
   
   
   maxDensity = Max(density)
   scaledDensity = BytScl(density, Min=0, Max=maxDensity)
   
   
   ; Set the number of contour levels
   if ask4 eq 1 then begin
    nlevels = (ceil(max(density)) - floor(min(density)))*2
    print, 'nlevels = ', nlevels
   endif else begin
    nlevels = (ceil(max(density)) - floor(min(density)))*100
    print, 'nlevels = ', nlevels
   endelse
   
   while nlevels gt 1000 do begin
    nlevels = nlevels/2
    print, 'nlevels =', nlevels
   endwhile
   

   
   ; Ascertain the step size needed to cover the data
   step = (max(density) - min(density))/nlevels
   
   ; Set the actual levels to be contoured
   levelvector = FINDGEN(nlevels)*step + min(density)
   
   lonarray = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0]
   biglonarray = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0, 195.0, 210.0, 225.0, 240.0, 255.0, 270.0, 285.0, 300.0, 315.0, 330.0, 345.0,360.0]
   if ask3 eq 2 then begin
     latarray = indgen((latmaxloc-latminloc), START = minx, INCREMENT = 5) - .5
   endif

   ct = COLORTABLE(74, /reverse)
   if ask3 eq 3 then begin 
     contourgitm = CONTOUR(density, Lons, Alts, $
              /FILL, RGB_TABLE = ct, YRANGE = yrange1, $
              XRANGE = xrange1, $
              YTITLE = 'Altitude (km)', $
              XTITLE = 'Longitude ($\circ$)', $
              TITLE = varname + ' Slice Average from ' + strtrim(lowestlatitude, 1) + ' to ' + strtrim(highestlatitude, 1) + ' degrees Latitude', $
              DIM = [1920, 1080], $
              C_VALUE = levelvector, $ 
              NAME = 'M-GITM',$
              XTICKVALUES = lonarray, $
              FONT_SIZE = 12, $
              XMINOR = 6, $
              MARGIN = [0.15, 0.15, 0.20, 0.15])
    endif
    
    if ask3 eq 2 then begin
      contourgitm = CONTOUR(density, Lats, Alts, $
        /FILL, RGB_TABLE = ct, YRANGE = yrange1, $
        XRANGE = xrange1, $
        YTITLE = 'Altitude (km)', $
        XTITLE = 'Latitude ($\circ$)', $
        TITLE = varname + ' Slice Average from ' + strtrim(lowestlongitude, 1) + ' to ' + strtrim(highestlongitude, 1) + ' degrees Longitude', $
        DIM = [1920, 1080], $
        C_VALUE = levelvector, $
        NAME = 'M-GITM',$
        XTICKVALUES = latarray, $
        FONT_SIZE = 12, $
        XTICKFONT_SIZE = 5.5, $
        XMINOR = 6, $
        MARGIN = [0.15, 0.15, 0.20, 0.15])
    endif
            
    local_time = ((biglonarray - meta.longsubsol)/15.0) + 12.0
    
    for g = 0,(n_elements(local_time)-1) do begin
      if local_time[g] gt 24 then begin
        local_time[g] = (local_time[g] - 24)
      endif
    endfor

    ;print, local_time
    if ask3 eq 3 then begin
      a_time = axis(0, $
        ;tickvalues = 0, $
        ticklen = (1.0/24.0), $
        tickname = strtrim(uint(local_time),1), $
        location=[0,min(yrange1)-25,0], $
        title='Local Time')
    endif
    
                    

    maxdenseint = floor(max(density))
    mindenseint = ceil(min(density))
    originaldenseintdifference = (maxdenseint - mindenseint)
  
    RF = 5.0
    
    
        
    tickv = []
    tickn = []
    tickstep = mindenseint
    
    if originaldenseintdifference ge 10 then begin
      maxdenseround = ceil(maxdenseint/RF)*RF
      originaldenseintdifference = ceil(originaldenseintdifference/RF)*RF
      tickstepdiff = ceil(ceil(originaldenseintdifference/10)/RF)*RF
      tickstep = (ceil(tickstep/RF)*RF)
      denseintdifference = tickstepdiff*1
    endif else begin
      maxdenseround = ceil(maxdenseint)
      denseintdifference = 1*originaldenseintdifference
      tickstepdiff = 1
    endelse
    
    while tickstep le maxdenseround do begin
      tickv = [tickv, tickstep]
      tickn = [tickn, strtrim(tickstep,1)]
      tickstep = tickstep + tickstepdiff
    endwhile
    
    
    im = image(density, RGB_table = ct, IMAGE_DIMENSIONS=[xcoordsteps,ycoordsteps], /overplot)
    c = COLORBAR(TARGET=im, TAPER = 0, $
            ORIENTATION = 1, $
            POSITION = [0.900, 0.25, 0.92, 0.75], $
            FONT_SIZE = 6, TICKVALUES = tickv, TICKNAME = tickn)
            
    ;;Single altitude Plots
    alt1 = 171.25
    alt2 = 201.25
    alt3 = 241.25
    
    altitude1 = density[*,where(alts eq (171.25))]
    altitude2 = density[*,where(alts eq (201.25))]
    altitude3 = density[*,where(alts eq (241.25))]
    

    if ask3 eq 3 then begin
      Pressure1plot = plot(Lons,altitude1, thick = 2, XRANGE = xrange1, XTICKVALUES = lonarray, YTITLE = 'Log Number Density (Log(#/m3))', $
        XTITLE = 'Longitude ($\circ$)', NAME = ('Altitude = ' + strtrim(alt1,1)), COLOR = 'Blue', DIM = [1920, 1080], $
        title = ('Number Density vs Longitude for ' + strtrim(variables[(ask4-1)],1) + ' at altitude = ' + strtrim(alt1,1) + ', ' + strtrim(alt2,1) + ', ' + strtrim(alt3,1) ))

      Pressure2plot = plot(Lons,altitude2, thick = 2, XRANGE = xrange1, XTICKVALUES = lonarray,  NAME = ('Altitude = ' + strtrim(alt2,1)), COLOR = 'Green', $
        /overplot)

      Pressure3plot = plot(Lons,altitude3, thick = 2, XRANGE = xrange1, XTICKVALUES = lonarray, NAME = ('Altitude = ' + strtrim(alt3,1)), COLOR = 'Gold', $
        /overplot)

      b_time = axis('X', ticklen = (1.0/24.0), tickname = strtrim(uint(local_time),1), location=[(min(altitude3)-1.2)], title='Local Time')

      legend = LEGEND(TARGET=[Pressure1plot,Pressure2plot,Pressure3plot], POSITION=[.9,.8], /AUTO_TEXT_COLOR)
    endif
    
    if ask3 eq 2 then begin
      Pressure1plot = plot(Lats,altitude1, thick = 2, XRANGE = xrange1, XTICKVALUES = latarray, YTITLE = 'Log Number Density (Log(#/m3))', $
        XTITLE = 'Latitude ($\circ$)', NAME = ('Altitude = ' + strtrim(alt1,1)), COLOR = 'Blue', DIM = [1920, 1080], $
        title = ('Number Density vs Latitude for ' + strtrim(variables[(ask4-1)],1) + ' at altitude = ' + strtrim(alt1,1) + ', ' + strtrim(alt2,1) + ', ' + strtrim(alt3,1) ), XTICKFONT_SIZE = 5.5)

      Pressure2plot = plot(Lats,altitude2, thick = 2, XRANGE = xrange1, XTICKVALUES = latarray,  NAME = ('Altitude = ' + strtrim(alt2,1)), COLOR = 'Green', $
        /overplot)

      Pressure3plot = plot(lats,altitude3, thick = 2, XRANGE = xrange1, XTICKVALUES = latarray, NAME = ('Altitude = ' + strtrim(alt3,1)), COLOR = 'Gold', $
        /overplot)

      legend = LEGEND(TARGET=[Pressure1plot,Pressure2plot,Pressure3plot], POSITION=[.9,.8], /AUTO_TEXT_COLOR)
    endif
    
;    ct = COLORTABLE(0, /reverse)
;    contourgitm2 = CONTOUR(RatioEqual1, Lons, Alts, $
;            /FILL, RGB_TABLE = ct, YRANGE = yrange1, $
;            XRANGE = xrange1, $
;            YTITLE = 'Altitude (km)', $
;            XTITLE = 'Longitude ($\circ$)', $
;            TITLE = varname + ' Slice Average from ' + strtrim(lowestlatitude, 1) + ' to ' + strtrim(highestlatitude, 1) + ' degrees Latitude', $
;            DIM = [1280, 720], $
;            C_VALUE = levelvector, $
;            NAME = 'M-GITM',$
;            XTICKVALUES = lonarray, $
;            FONT_SIZE = 6, $
;            XMINOR = 6, $
;            MARGIN = [0.15, 0.15, 0.20, 0.15], /overplot)
            
    CD, picsavloc
    Contourgitm.save, (Variables[(ask4-1)] + ' contour plot altitudes ' + strtrim(miny,1) + '-' + strtrim(maxy,1) + ' latitudes ' + strtrim(minz,1) + '-' + strtrim(maxz,1) + '(altitude scale).png')
    pressure1plot.save, (Variables[(ask4-1)] + ' at constant altitudes ' + strtrim(170,1) + ', ' + strtrim(200,1) + ', ' + strtrim(240,1) +  '.png')
    Processedfilename = outputsavloc + name + ' Altitude ' + Variables[ask4-1] + ' at ' + strtrim(namez,1) + ' ' + strtrim(FIX(minz),1) + '-' + strtrim(FIX(maxz),1) + ' ' + strtrim(namex,1) + ' ' + strtrim(Fix(minx),1) + '-' + strtrim(FIX(maxx),1) + '.sav'
    if name ne 'n' then begin
      save, ask1, ask2, ask3, ask4, varname, latarray, lonarray,  local_time, small_alts, small_p_interp, lons, lats, xcoordsteps, ycoordsteps, biglonarray, lowestlatitude, highestlatitude, lowestlongitude, highestlongitude, density, xrange1, yrange1, p_interp, filename = Processedfilename
      print, 'saved'
    endif
    save, density, filename = Processedfilename
stop        
END