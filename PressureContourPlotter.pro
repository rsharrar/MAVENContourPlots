PRO DataCubeContourPlotterPressure
;;Save pictures location
picsavloc = "C:\Users\rshar\Pictures"
;;Save output location
outputsavloc = 'C:\Users\rshar\IDLWorkspace83\Default\Current Data\April 2017\2D Plot Test Area\Processed Data\'
;;Pressure Version

file = dialog_pickfile(path = 'C:\Users\rshar\IDLWorkspace83\Default\Current Data\April 2017\2D Plot Test Area')
restore, file

Variables = ['Temperature','Co2', 'Co','N2','O','Co-Co2 Ratio', 'N2-Co2 Ratio','O-Co2 Ratio']
VarNameLength = n_elements(Variables)

print, '1: Pressure'
print, '2: Longitude'
print, '3: Latitude'

ask1 = float(ask('the number corresponding to the variable that you would like to be the Y coordinate '))
ask2 = float(ask('the number corresponding to the variable that you would like to be the X coordinate '))
ask3 = float(ask('the number corresponding to the variable that you would like to be the Z coordinate '))

for v = 0,(VarNameLength - 1) do begin
  print, strtrim(v+1,1), ': ', Variables[v]
endfor

ask4 = float(ask('the corresponding number for the Variable you would like to plot'))

miny = float(ask('the minimum altitude you would like to look at'))
maxy = float(ask('the maximum altitude you would like to look at'))


if ask3 eq 3 then begin
  minz = float(ask('the minimum Latitude in your range'))
  maxz = float(ask('the maximum Latitude in your range'))
  namez = 'Lat'
  namex = 'Lon
endif

if ask3 eq 2 then begin
  minz = float(ask('the minimum Longitude in your range'))
  maxz = float(ask('the maximum Longitude in your range'))
  namez = 'Lon'
  namex = 'Lat'
endif

name = string(ask('a name for this run if you would like to save the processed data. Otherwise, enter "n" '))

;ask1 = 1.0
;ask2 = 2.0
;ask3 = 3.0
;ask4 = 1.0
;miny = 103
;maxy = 250.0
;minz = -15
;maxz = 15

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

;ycoordsteps = (maxy - miny)/binalt + 1
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
  lowestlongitude = 0
  highestlongitude = 0
endif

if ask3 eq 2 then begin
  lowestlongitude = minz
  Zdim = lowestlongitude*1
  Zdimp = 1*Zdim
  highestlongitude = maxz
  lowestlatitude = 0
  highestlatitude = 0
endif

;;Finding common pressure range for latitude slices
;;Setup the Y Coordinate according to user input
if ask1 eq 1 then begin
  ;yrange1 = [minalt,maxalt]
  ybinsize = binalt
  altminloc = fix((miny - minalt)/binalt)
  altmaxloc = fix((maxy - minalt)/binalt)
endif

if ask1 eq 2 then begin
  ;yrange1 = [minlon,maxlon]
  ybinsize = binlon
  lonminloc = fix((miny - minlon)/binlon)
  lonmaxloc = fix((maxy - minlon)/binlon)
endif

if ask1 eq 3 then begin
  ;yrange1 = [minlat,maxlat]
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

minpressures = []
maxpressures = []

minpressureloc = []
maxpressureloc = []

;;Pressure Interp for Latitude
if ask3 eq 3 then begin
  latitude_slices = abs((lowestlatitude - highestlatitude)/binlat)
  
  for q = 0,(latitude_slices) do begin
  ;;Setup the Z coordinate according to user input
    ;Debuggin code: print, 'Interpolating slice at', Zdim
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
  
    minpressure_i = min(ndensitys.P[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    maxpressure_i = max(ndensitys.P[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
  
    minpressures = [[minpressures],[minpressure_i]]
    maxpressures = [[maxpressures],[maxpressure_i]]
  
    for a = lonminloc,lonmaxloc do begin
      for b = (altminloc-39),(altmaxloc-39) do begin
        if (ndensitys.P[a,latminloc:latmaxloc,b] eq minpressure_i) then begin
          minpressureloc = [[minpressureloc],[b]]
        endif
        
        if (ndensitys.P[a,latminloc:latmaxloc,b] eq maxpressure_i) then begin
          maxpressureloc = [[maxpressureloc],[b]]
        endif
        
      endfor
    endfor
    
    Zdim = Zdim + binlat
  
  endfor
endif

;;Pressure Interpol for Longitude
if ask3 eq 2 then begin
  longitude_slices = abs((lowestlongitude - highestlongitude)/binlat)

  for q = 0,(longitude_slices) do begin
    ;;Setup the Z coordinate according to user input
    ;Debuggin code: print, 'Interpolating slice at', Zdim
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

    minpressure_i = min(ndensitys.P[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])
    maxpressure_i = max(ndensitys.P[lonminloc:lonmaxloc,latminloc:latmaxloc,(altminloc-39):(altmaxloc-39)])

    minpressures = [[minpressures],[minpressure_i]]
    maxpressures = [[maxpressures],[maxpressure_i]]

    for a = latminloc,latmaxloc do begin
      for b = (altminloc-39),(altmaxloc-39) do begin
        if (ndensitys.P[lonminloc:lonmaxloc,a,b] eq minpressure_i) then begin
          minpressureloc = [[minpressureloc],[b]]
        endif

        if (ndensitys.P[lonminloc:lonmaxloc,a,b] eq maxpressure_i) then begin
          maxpressureloc = [[maxpressureloc],[b]]
        endif

      endfor
    endfor

    Zdim = Zdim + binlon

  endfor
endif

alts = meta.altitude[(min(maxpressureloc)+39):(max(minpressureloc)+39)]

;;Making Interpolation Vector
minpressure = max(minpressures)
maxpressure = min(maxpressures)

logminp = log(minpressure)
logmaxp = log(maxpressure)

points     = abs((floor(logmaxp)-ceil(logminp))*100)
p_interp = findgen(points+1)/100+ceil(logminp)
p_interp = reverse(p_interp)
ycoordsteps = float(N_ELEMENTS(p_interp))
yrange1 = [max(p_interp),min(p_interp)]

;;Setting up for loop
X = []
factor = 0
top = MAKE_ARRAY(xcoordsteps, ycoordsteps, /float, VALUE = 0.0)
top2 = MAKE_ARRAY(xcoordsteps, ycoordsteps, /float, VALUE = 0.0)
zdim = 1*zdimp

alt_interp_array = []

if ask3 eq 3 then begin
  
  for n = 0,(latitude_slices) do begin
    densityslice_i = []
    densityslice2_i = []
    ;Debuggin code: print, 'Averaging slice at', Zdim
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
      for i = (lonminloc),(lonmaxloc) do begin
        minislice = abs(interpol(reform(temperature.Tn[i,latminloc,(altminloc-39):(altmaxloc-39)]),log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp)))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'Temperature'
    endif
  
    if ask4 eq 2 then begin
      for i = (lonminloc),(lonmaxloc) do begin
        minislice = abs(interpol(log(reform(ndensitys.Co2[i,latminloc,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp)))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'Co2'
    endif
    
    if ask4 eq 3 then begin
      for i = (lonminloc),(lonmaxloc) do begin
        minislice = interpol(log(reform(ndensitys.Co[i,latminloc,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'Co'
    endif
    
    if ask4 eq 4 then begin
      for i = (lonminloc),(lonmaxloc) do begin
        minislice = interpol(log(reform(ndensitys.N2[i,latminloc,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'N2'
    endif
    
    if ask4 eq 5 then begin
      for i = (lonminloc),(lonmaxloc) do begin
        minislice = interpol(log(reform(ndensitys.O[i,latminloc,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'Oxygen'
    endif
    
    if ((ask4 ge 6) and (ask4 le 8)) then begin
      varname = Variables[(ask4-1)]
      for i = (lonminloc),(lonmaxloc) do begin
        minislice = abs(interpol(log(reform(ndensitys.Co2[i,latminloc,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp)))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      if ask4 eq 6 then begin
        for i = (lonminloc),(lonmaxloc) do begin
          minislice = abs(interpol(log(reform(ndensitys.Co[i,latminloc,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp)))
          densityslice2_i = [[densityslice2_i],[minislice]]
        endfor
      endif
      if ask4 eq 7 then begin
        for i = (lonminloc),(lonmaxloc) do begin
          minislice = abs(interpol(log(reform(ndensitys.N2[i,latminloc,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp)))
          densityslice2_i = [[densityslice2_i],[minislice]]
        endfor
      endif
      if ask4 eq 8 then begin
        for i = (lonminloc),(lonmaxloc) do begin
          minislice = abs(interpol(log(reform(ndensitys.O[i,latminloc,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp)))
          densityslice2_i = [[densityslice2_i],[minislice]]
        endfor
      endif
    endif
    
    factor_i = cos(Zdim*!dtor)
    
    densityslice_i = 1*transpose(densityslice_i)
    
    if (ask4 ge 6) and (ask4 le 8) then begin
      densityslice2_i = 1*transpose(densityslice2_i)
    endif
    
    top_i = factor_i*densityslice_i
    
    if (ask4 ge 6) and (ask4 le 8) then begin
      top2_i = factor_i*densityslice2_i
    endif
    
    factor = factor + factor_i
    top = top + top_i
    if (ask4 ge 6) and (ask4 le 8) then begin
      top2 = top2 +top2_i
    endif
    
    Zdim = Zdim + binlat
    
    for i = (lonminloc),(lonmaxloc) do begin
      alt_interpolated2_i = interpol(alts,log(reform(ndensitys.P[i,latminloc,(altminloc-39):(altmaxloc-39)])),(p_interp))
      alt_interp_array = [[alt_interp_array],[alt_interpolated2_i]]
      if (max(alt_interpolated2_i) le maxy) and (min(alt_interpolated2_i) ge miny) then begin
        alt_interpolated2 = alt_interpolated2_i
      endif
    endfor
    
  endfor
endif

;;Average over Longitude Slices
if ask3 eq 2 then begin
  for n = 0,(longitude_slices) do begin
    densityslice_i = []
    densityslice2_i = []
    ;Debuggin code: print, 'Averaging slice at', Zdim
    ;;Setup the Z coordinate according to user input
    lonminloc = fix(Zdim - minlon)/binlon
    lonmaxloc = fix(Zdim - minlon)/binlon
    


    ; Set up variables for the plot. For now, using Altitude as the Y coordinate and Longitude as the X coordinate
    if ask4 eq 1 then begin
      for i = (latminloc),(latmaxloc) do begin
        minislice = abs(interpol(reform(temperature.Tn[lonminloc,i,(altminloc-39):(altmaxloc-39)]),log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp)))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'Temperature'
    endif

    if ask4 eq 2 then begin
      for i = (latminloc),(latmaxloc) do begin
        minislice = abs(interpol(log(reform(ndensitys.co2[lonminloc,i,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp)))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'Co2'
    endif

    if ask4 eq 3 then begin
      for i = (latminloc),(latmaxloc) do begin
        minislice = abs(interpol(log(reform(ndensitys.co[lonminloc,i,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp)))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'Co'
    endif

    if ask4 eq 4 then begin
      for i = (latminloc),(latmaxloc) do begin
        minislice = abs(interpol(log(reform(ndensitys.n2[lonminloc,i,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp)))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'N2'
    endif

    if ask4 eq 5 then begin
      for i = (latminloc),(latmaxloc) do begin
        minislice = abs(interpol(log(reform(ndensitys.O[lonminloc,i,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp)))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      varname = 'Oxygen'
    endif

    if ((ask4 ge 6) and (ask4 le 8)) then begin
      varname = Variables[(ask4-1)]
      for i = (latminloc),(latmaxloc) do begin
        minislice = abs(interpol(log(reform(ndensitys.co2[lonminloc,i,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp)))
        densityslice_i = [[densityslice_i],[minislice]]
      endfor
      if ask4 eq 6 then begin
        for i = (latminloc),(latmaxloc) do begin
          minislice = abs(interpol(log(reform(ndensitys.co[lonminloc,i,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp)))
          densityslice2_i = [[densityslice2_i],[minislice]]
        endfor
      endif
      if ask4 eq 7 then begin
        for i = (latminloc),(latmaxloc) do begin
          minislice = abs(interpol(log(reform(ndensitys.n2[lonminloc,i,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp)))
          densityslice2_i = [[densityslice2_i],[minislice]]
        endfor
      endif
      if ask4 eq 8 then begin
        for i = (latminloc),(latmaxloc) do begin
          minislice = abs(interpol(log(reform(ndensitys.O[lonminloc,i,(altminloc-39):(altmaxloc-39)])),log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp)))
          densityslice2_i = [[densityslice2_i],[minislice]]
        endfor
      endif
    endif

    factor_i = 1

    densityslice_i = 1*transpose(densityslice_i)

    if (ask4 ge 6) and (ask4 le 8) then begin
      densityslice2_i = 1*transpose(densityslice2_i)
    endif

    top_i = factor_i*densityslice_i

    if (ask4 ge 6) and (ask4 le 8) then begin
      top2_i = factor_i*densityslice2_i
    endif

    factor = factor + factor_i
    top = top + top_i
    if (ask4 ge 6) and (ask4 le 8) then begin
      top2 = top2 +top2_i
    endif

    Zdim = Zdim + binlon

    for i = (latminloc),(latmaxloc) do begin
      alt_interpolated2_i = interpol(alts,log(reform(ndensitys.P[lonminloc,i,(altminloc-39):(altmaxloc-39)])),(p_interp))
      alt_interp_array = [[alt_interp_array],[alt_interpolated2_i]]
      if ((max(alt_interpolated2_i) le maxy) and (max(alt_interpolated2_i) ge (maxy-20))) and ((min(alt_interpolated2_i) ge miny) and (min(alt_interpolated2_i) le miny + 20)) then begin
        alt_interpolated2 = alt_interpolated2_i
      endif
    endfor

  endfor
endif
print, alt_interpolated2

  dens1 = top/factor
  dens2 = top2/factor
  
  
  if ((ask4 ge 6) and (ask4 le 8)) then begin
   X = dens2 - dens1
   print, 'Ratio'
  endif else begin
   X = dens1
   print, 'No Ratio'
  endelse
   
   lons = meta.longitude[lonminloc:lonmaxloc]
   minlonaxis = min(lons)
   maxlonaxis = max(lons)
   
   lats = meta.latitude[latminloc:latmaxloc]
   minlataxis = min(lats)
   maxlataxis = max(lats)
   
   if ask3 eq 3 then begin
     xrange1 = [minlonaxis,maxlonaxis]
   endif
   
   if ask3 eq 2 then begin
     xrange1 = [minlataxis,maxlataxis]
   endif
   
   ; Create the density plot by binning the data into a 2D histogram.
   density = X
   
   
   
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
   endif else begin
    latarray = [0]
   endelse
   

   ct = COLORTABLE(74, /reverse)
   if ask3 eq 3 then begin
      contourgitm = CONTOUR(density, Lons, p_interp, $
            /FILL, RGB_TABLE = ct, YRANGE = yrange1, $
            XRANGE = xrange1, $
            YTITLE = 'Log Pressure (Log(Pascals))', $
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
      contourgitm = CONTOUR(density, Lats, p_interp, $
        /FILL, RGB_TABLE = ct, YRANGE = yrange1, $
        XRANGE = xrange1, $
        YTITLE = 'Log Pressure (Log(Pascals))', $
        XTITLE = 'Latitude ($\circ$)', $
        TITLE = varname + ' Slice Average from ' + strtrim(lowestlongitude, 1) + ' to ' + strtrim(highestlongitude, 1) + ' degrees Longitude', $
        DIM = [1920, 1080], $
        C_VALUE = levelvector, $
        NAME = 'M-GITM',$
        XTICKVALUES = latarray, $
        FONT_SIZE = 12, $
        XMINOR = 6, $
        XTICKFONT_SIZE = 5.5, $
        MARGIN = [0.15, 0.15, 0.20, 0.15])
    endif
            
    if ask3 eq 3 then begin
      local_time = ((biglonarray - meta.longsubsol)/15.0) + 12.0
    endif
    
    for g = 0,(n_elements(local_time)-1) do begin
      if local_time[g] gt 24 then begin
        local_time[g] = (local_time[g] - 24)
      endif
    endfor

    ;print, local_time
    
    small_alts = reverse([alt_interpolated2[(points-375)],alt_interpolated2[(points-350)],alt_interpolated2[(points-325)],alt_interpolated2[(points-300)],alt_interpolated2[(points-275)],alt_interpolated2[(points-250)],alt_interpolated2[(points-225)],alt_interpolated2[(points-200)],alt_interpolated2[(points-175)],alt_interpolated2[(points-150)],alt_interpolated2[(points-125)],alt_interpolated2[(points-100)],alt_interpolated2[(points-75)],alt_interpolated2[(points-50)],alt_interpolated2[(points-25)],alt_interpolated2[(points-0)]])
    small_p_interp = [p_interp[(points-375)],p_interp[(points-350)],p_interp[(points-325)],p_interp[(points-300)],p_interp[(points-275)],p_interp[(points-250)],p_interp[(points-225)],p_interp[(points-200)],p_interp[(points-175)],p_interp[(points-150)],p_interp[(points-125)],p_interp[(points-100)],p_interp[(points-75)],p_interp[(points-50)],p_interp[(points-25)],p_interp[(points-0)]]
    
    if ask3 eq 3 then begin
      ;alt_interpol_axis = axis('Y', tickvalues = small_p_interp, tickname = strtrim(uint(small_alts),1), location=[(max(lons))], title='Altitude', textpos = 1)
      a_time = axis('X', ticklen = (1.0/24.0), tickname = strtrim(uint(local_time),1), location=[(min(p_interp))+min(p_interp)*.075], title='Local Time')
    endif
    
    if ask3 eq 2 then begin
      ;alt_interpol_axis = axis('Y', tickvalues = small_p_interp, tickname = strtrim(uint(small_alts),1), location=[(max(lats)+3)], title='Altitude', textpos = 1, TICKFONT_SIZE = 6)
    endif
    
    rangemin = (round(min(density)*4))/4
    rangemax = (round(max(density)*4))/4
    
    if ((ask4 lt 6)) then begin
      range1 = [rangemin,rangemax]
    endif else begin
      range1 = [min(density),max(density)]  
    endelse
    
    c = COLORBAR(RGB_table = ct, TAPER = 0, $
      ORIENTATION = 1, $
      POSITION = [0.900, 0.25, 0.92, 0.75], $
      FONT_SIZE = 6, RANGE = range1)
    
    if ((ask4 ge 6) and (ask4 le 8)) then begin
          RatioEqual1 = MAKE_ARRAY(xcoordsteps, ycoordsteps, /float, VALUE = log(-1))
          RatioEqual1[WHERE(((density gt .98) and (density lt 1.02)), /NULL)] = 1
         endif
         
    ;;Single Pressure Plots
    pressure1 = density[*,where(p_interp eq (ceil(logminp) + 0))]
    pressure2 = density[*,where(p_interp eq (ceil(logminp) + 1))]
    pressure3 = density[*,where(p_interp eq (ceil(logminp) + 2))]
    pressure4 = density[*,where(p_interp eq (ceil(logminp) + 3))]
    
    if ask3 eq 3 then begin
      Pressure1plot = plot(Lons,pressure1, thick = 2, XRANGE = xrange1, XTICKVALUES = lonarray, YTITLE = 'Log Number Density (Log(#/m3))', $
        XTITLE = 'Longitude ($\circ$)', NAME = ('Pressure = ' + strtrim(ceil(logminp) + 0,1)), COLOR = 'Blue', DIM = [1920, 1080], $
        title = ('Number Density vs Longitude for ' + strtrim(variables[(ask4-1)],1) + ' at log(pressure) = ' + strtrim(ceil(logminp) + 0,1) + ', ' + strtrim(ceil(logminp) + 1,1) + ', ' + strtrim(ceil(logminp) + 2,1) + ', ' + strtrim(ceil(logminp) + 3,1)))

      Pressure2plot = plot(Lons,pressure2, thick = 2, XRANGE = xrange1, XTICKVALUES = lonarray,  NAME = ('Pressure = ' + strtrim(ceil(logminp) + 1,1)), COLOR = 'Green', $
        /overplot)

      Pressure3plot = plot(Lons,pressure3, thick = 2, XRANGE = xrange1, XTICKVALUES = lonarray, NAME = ('Pressure = ' + strtrim(ceil(logminp) + 2,1)), COLOR = 'Gold', $
        /overplot)

      Pressure4plot = plot(Lons,pressure4, thick = 2, XRANGE = xrange1, XTICKVALUES = lonarray, NAME = ('Pressure = ' + strtrim(ceil(logminp) + 3,1)), COLOR = 'Red', $
        /overplot)

      b_time = axis('X', ticklen = (1.0/24.0), tickname = strtrim(uint(local_time),1), location=[(min(pressure4))*.835], title='Local Time')

      legend = LEGEND(TARGET=[Pressure1plot,Pressure2plot,Pressure3plot,Pressure4plot], POSITION=[.9,.8], /AUTO_TEXT_COLOR)
    endif
    
    if ask3 eq 2 then begin
      Pressure1plot = plot(Lats,pressure1, thick = 2, XRANGE = xrange1, XTICKVALUES = latarray, YTITLE = 'Log Number Density (Log(#/m3))', $
        XTITLE = 'Latitude ($\circ$)', NAME = ('Pressure = ' + strtrim(ceil(logminp) + 0,1)), COLOR = 'Blue', DIM = [1920, 1080], $
        title = ('Number Density vs Latitude for ' + strtrim(variables[(ask4-1)],1) + ' at log(pressure) = ' + strtrim(ceil(logminp) + 0,1) + ', ' + strtrim(ceil(logminp) + 1,1) + ', ' + strtrim(ceil(logminp) + 2,1) + ', ' + strtrim(ceil(logminp) + 3,1)), XTICKFONT_SIZE = 5.5)

      Pressure2plot = plot(Lats,pressure2, thick = 2, XRANGE = xrange1, XTICKVALUES = latarray,  NAME = ('Pressure = ' + strtrim(ceil(logminp) + 1,1)), COLOR = 'Green', XTICKFONT_SIZE = 5.5, $
        /overplot)

      Pressure3plot = plot(Lats,pressure3, thick = 2, XRANGE = xrange1, XTICKVALUES = latarray, NAME = ('Pressure = ' + strtrim(ceil(logminp) + 2,1)), COLOR = 'Gold', XTICKFONT_SIZE = 5.5, $
        /overplot)

      Pressure4plot = plot(Lats,pressure4, thick = 2, XRANGE = xrange1, XTICKVALUES = latarray, NAME = ('Pressure = ' + strtrim(ceil(logminp) + 3,1)), COLOR = 'Red', XTICKFONT_SIZE = 5.5, $
        /overplot)

      legend = LEGEND(TARGET=[Pressure1plot,Pressure2plot,Pressure3plot,Pressure4plot], POSITION=[.9,.8], /AUTO_TEXT_COLOR)
    endif
     
;   ct2 = COLORTABLE(0, /reverse)
;       contourgitm2 = CONTOUR(RatioEqual1, Lons, p_interp, $
;               /FILL, RGB_TABLE = ct2, YRANGE = yrange1, $
;               XRANGE = xrange1, $
;               YTITLE = 'Log Pressure (Log(Pascals))', $
;               XTITLE = 'Longitude ($\circ$)', $
;               TITLE = varname + ' Slice Average from ' + strtrim(lowestlatitude, 1) + ' to ' + strtrim(highestlatitude, 1) + ' degrees Latitude', $
;               DIM = [1280, 720], $
;               C_VALUE = levelvector, $
;               NAME = 'M-GITM',$
;               XTICKVALUES = lonarray, $
;               FONT_SIZE = 6, $
;               XMINOR = 6, $
;               MARGIN = [0.15, 0.15, 0.20, 0.15], /overplot)
;       
    
          
    
            
    CD, picsavloc
    Contourgitm.save, (Variables[(ask4-1)] + ' contour plot altitudes ' + strtrim(miny,1) + '-' + strtrim(maxy,1) + ' latitudes ' + strtrim(minz,1) + '-' + strtrim(maxz,1) + '.png')
    pressure1plot.save, (Variables[(ask4-1)] + '  at constant pressures ' + strtrim(ceil(logminp) + 0,1) + ', ' + strtrim(ceil(logminp) + 1,1) + ', ' + strtrim(ceil(logminp) + 2,1) + ', ' + strtrim(ceil(logminp) + 3,1) + '.png')
    Processedfilename = outputsavloc + name + ' Pressure Scale ' + Variables[ask4-1] + ' at ' + strtrim(namez,1) + ' ' + strtrim(FIX(minz),1) + '-' + strtrim(FIX(maxz),1) + ' ' + strtrim(namex,1) + ' ' + strtrim(Fix(minx),1) + '-' + strtrim(FIX(maxx),1) + '.sav'
    if name ne 'n' then begin
       save, ask1, ask2, ask3, ask4, varname, latarray, lonarray,  local_time, small_alts, small_p_interp, lons, lats, xcoordsteps, ycoordsteps, biglonarray, lowestlatitude, highestlatitude, lowestlongitude, highestlongitude, density, xrange1, yrange1, p_interp, filename = Processedfilename  
       print, 'saved'
    endif
      
print, 'Done'
END