pro datacubeprocesseddataaverager
  file = dialog_pickfile(path = 'C:\Users\rshar\IDLWorkspace83\Default\Current Data\April 2017\2D Plot Test Area\Processed Data', /MULTIPLE_FILES)
  
  numfiles = n_elements(file)
  restore, file[0]
  totaldensity = density
  
for n = 1,(numfiles-1) do begin
  restore, file[n]
  totaldensity += density
endfor

density = (totaldensity/numfiles)

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

for g = 0,(n_elements(local_time)-1) do begin
  if local_time[g] gt 24 then begin
    local_time[g] = (local_time[g] - 24)
  endif
endfor

if ask3 eq 3 then begin
  ;alt_interpol_axis = axis('Y', tickvalues = small_p_interp, tickname = strtrim(uint(small_alts),1), location=[(max(lons)+3)], title='Altitude', textpos = 1)
  a_time = axis('X', ticklen = (1.0/24.0), tickname = strtrim(uint(local_time),1), location=[(min(p_interp)-.5)], title='Local Time')
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
  FONT_SIZE = 12, RANGE = range1)

if ((ask4 ge 6) and (ask4 le 8)) then begin
  RatioEqual1 = MAKE_ARRAY(xcoordsteps, ycoordsteps, /float, VALUE = log(-1))
  RatioEqual1[WHERE(((density gt .98) and (density lt 1.02)), /NULL)] = 1
endif  
stop
end