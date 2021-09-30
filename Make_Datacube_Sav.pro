filelist = file_search('C:\Users\rshar\IDLWorkspace83\Default\Current Data\April 2017\2D Plot Test Area\Bin\3DALL*.bin')

if n_elements(ls) eq 0 then ls = 0
ls = fix(ask('LS: ',tostr(ls)))

if n_elements(subsolarlon) eq 0 then subsolarlon = 0
subsolarlon = fix(ask('subsolar longitude: ',tostr(subsolarlon)))

nFiles = n_elements(filelist)
if nFiles eq 1 then begin
	filename = filelist
endif else begin
	display, filelist
	if n_elements(fn) eq 0 then fn = 0
	fn = fix(ask('which file', tostr(fn)))
	filename = filelist(fn)
endelse	

read_thermosphere_file, filename, nvars, nalts, nlats, nlons, $
      vars, data, rb, cb, bl_cnt,time,version
      

latitude = data(1,0,2:nlats-3,0)

lon = reform(data(0,2:nlons-3,0,0))
lat = reform(data(1,0,2:nlats-3,0))
alt = reform(data(2,0,0,2:nalts-3))/1000.

gitm_alt_min = min(alt)
gitm_alt_max = max(alt)
;; initialize the user_alt_min variables
user_alt_min = 70.0
user_alt_max = 295.0

user_alt_min = double(fix(ask('Set Minimum Altitude (km)',tostr(user_alt_min))))
user_alt_max = double(fix(ask('Set Maximum Altitude (km)',tostr(user_alt_max))))

;;; Quick Error Checking to see if the user
;;; selected an altitude outside the model range
if (user_alt_min lt gitm_alt_min) then begin
    print, 'Warning: Minimum altitude selected below model limit, correcting.'
    user_alt_min = gitm_alt_min
endif 

if (user_alt_max gt gitm_alt_max) then begin
    print, 'Warning: Maximum altitude selected above mode limit, correcting.'
    user_alt_max = gitm_alt_max
endif 

;;; This sets the desired altitude range of the output
;;; so these are actually our aerodetic ranges
ialtmin = min(where(alt ge gitm_alt_min))
ialtmax = max(where(alt le gitm_alt_max))
naltsprint = ialtmax - ialtmin + 1

iSZA = where(vars eq 'SolarZenithAngle')
;
iCO2 = where(vars eq '[CO!D2!N]')
iCO = where(vars eq '[CO]')
iN2 = where(vars eq '[N!D2!N]')
iO2 = where(vars eq '[O!D2!N]')
iO = where(vars eq '[O]')
;
iOP = where(vars eq '[O!U+!N]')
iO2P = where(vars eq '[O!D2!U+!N]')
iCO2P = where(vars eq '[CO!D2!U+!N]')
ie = where(vars eq '[e-]')
;
iVeast = where(vars eq 'V!Dn!N(east)')
iVnorth = where(vars eq 'V!Dn!N(north)')
iVup = where(vars eq 'V!Dn!N(up)')
;
iTn = where(vars eq 'Temperature')
iTi = where(vars eq 'iTemperature')
iTe = where(vars eq 'eTemperature')
;
iQO = where(vars eq 'EUVIonizationRate(O!U+!N')
iQO2 = where(vars eq 'EUVIonizationRate(O!D2!U+!N')
iQCO2 = where(vars eq 'EUVIonizationRate(CO!D2!U+!N')
iQN2 = where(vars eq 'EUVIonizationRate(N!D2!U+!N')
iQNO = where(vars eq 'EUVIonizationRate(NO!U+!N')

nCO2 = reform(data(iCO2,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
nCO = reform(data(iCO,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
nN2 = reform(data(iN2,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
nO2 = reform(data(iO2,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
nO = reform(data(iO,2:nlons-3,2:nlats-3,ialtmin:ialtmax))

nOP = reform(data(iOP,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
nO2P = reform(data(iO2P,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
nCO2P = reform(data(iCO2P,2:nlons-3,2:nlats-3,ialtmin:ialtmax))

if ie ge 0 then begin
	n_e = reform(data(ie,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
endif else begin
	n_e = nOP + nO2P + nCO2P
endelse

vEast = reform(data(iVeast,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
vNorth = reform(data(iVNorth,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
vUp= reform(data(iVup,2:nlons-3,2:nlats-3,ialtmin:ialtmax))

Tn = reform(data(iTn,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
Ti = reform(data(iTi,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
Te = reform(data(iTe,2:nlons-3,2:nlats-3,ialtmin:ialtmax))

QOP1 = reform(data(iQO,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
QO2 = reform(data(iQO2,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
QN2 = reform(data(iQN2,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
QCO2 = reform(data(iQCO2,2:nlons-3,2:nlats-3,ialtmin:ialtmax))
QOP2 = reform(data(iQNO,2:nlons-3,2:nlats-3,ialtmin:ialtmax))

;;; JMB:  Note that some of these variables may be missing
if(iSZA lt 0) then begin
   SZA = fltarr(nlons-4,nlats-4)
   SZA(0:nlons-5,0:nlats-5) = 1.0e-02
endif else begin
   sza = reform(data(iSZA,2:nLons-3,2:nlats-3,ialtmin:ialtmax))
endelse

if(iQO lt 0) then begin
   QOP1 = fltarr(nlons-4,nlats-4,naltsprint)
   QOP1(0:nlons-5,0:nlats-5,0:naltsprint-1) = 1.0e-02
endif else begin
   QOP1 = reform(data(iQO,2:nLons-3,2:nlats-3,ialtmin:ialtmax))
endelse

if(iQO2 lt 0) then begin
   QO2 = fltarr(nlons-4,nlats-4,naltsprint)
   QO2(0:nlons-5,0:nlats-5,0:naltsprint-1) = 1.0e-02
endif else begin
   QO2 = reform(data(iQO2,2:nLons-3,2:nlats-3,ialtmin:ialtmax))
endelse

if(iQN2 lt 0) then begin
   QN2 = fltarr(nlons-4,nlats-4,naltsprint)
   QN2(0:nlons-5,0:nlats-5,0:naltsprint-1) = 1.0e-02
endif else begin
   QN2 = reform(data(iQN2,2:nLons-3,2:nlats-3,ialtmin:ialtmax))
endelse

if(iQCO2 lt 0) then begin
   QCO2 = fltarr(nlons-4,nlats-4,naltsprint)
   QCO2(0:nlons-5,0:nlats-5,0:naltsprint-1) = 1.0e-02
endif else begin
   QCO2 = reform(data(iQCO2,2:nLons-3,2:nlats-3,ialtmin:ialtmax))
endelse

iQNO = where(vars eq 'EUVIonizationRate(NO!U+!N')
if(iQNO lt 0) then begin
   QOP2 = fltarr(nlons-4,nlats-4,naltsprint)
   QOP2(0:nlons-5,0:nlats-5,0:naltsprint-1) = 1.0e-02
endif else begin
   QOP2 = reform(data(iQNO,2:nLons-3,2:nlats-3,ialtmin:ialtmax))
endelse

l = strpos(filename,'.bin',/reverse_search,/reverse_offset)-13
date = strmid(filename(0),l,13)
cyear = strmid(date,0,2)
if cyear gt 50 then year = 1900 + fix(cyear) else year = 2000 + fix(cyear)
month = fix(strmid(date,2,2))
day = fix(strmid(date,4,2))
hour = fix(strmid(date,7,2))
min = fix(strmid(date,9,2))
sec = fix(strmid(date,11,2))

;--------------------
; CONVERT THE ALTITUDE AND LATITUDE TO AEROCENTRIC:  JMB 06/22/2017
;
; Our current geometry requires
; (1) naltsprint (ialtmax - ialtmin + 1)
; (2) nlats - 4
; (3) nlons - 4
;
; The Aerocentric-Aerodetic coordinates
; cause convolution between latitude and altitude
; longitude is not involved and remains unchanged 
; in the transformation
; Hence, we set iLon = 0

AeroGITMAlts = fltarr(nlats-4,naltsprint)
AeroGITMLats = fltarr(nlats-4,naltsprint)

iLon = 0
for iAlt = ialtmin,ialtmax do begin
;  print, 'iAlt =', iAlt
  for iLat = 0,nlats-5 do begin
;     print, 'iLat =', iLat
     geocentric = [lat[iLat]*(180.0/!pi), lon[iLon]*(180.0/!pi), alt[iAlt]]
;     print, 'geocentric = ', geocentric
     newcoords = geo2geodetic(geocentric,PLANET='Mars')
     AeroGITMAlts[iLat,iAlt-ialtmin] = newcoords[2] 
     AeroGITMLats[iLat,iAlt-ialtmin] = newcoords[0] 

;     print, 'alts, AeroGITMAlts = ', alt[iAlt], newcoords[2]
;     print, 'lats, AeroGITMLats = ', lat[iLat]*(180.0/!pi), newcoords[0]
  endfor 
endfor 

; next, interpolate our aerodetic lats, assume that there is no
; variation in altitude, which is mostly correct
; SPLINE SYNTAX:
; Result = SPLINE(X,Y,T, /DOUBLE).
;     X = original indep. variable  (i.e., latitude)
;     Y = original dep. variable 
;     T = new indep. variable points (i.e., geodetic latitude)

; for our interpolation
; X = original Aerodetic Latitudes (assumed constant with altitude)
; Y = our chosen variable (e.g. Tn or Ti, etc.)
; T = our desired Aerodetic grid

; non-log variables
NewLats = fltarr(nlats-4)
NewAlts = fltarr(nlats-4,naltsprint)

NewTn = fltarr(nlons-4,nlats-4,naltsprint)
NewTe = fltarr(nlons-4,nlats-4,naltsprint)
NewTi = fltarr(nlons-4,nlats-4,naltsprint)
;
NewnCO2 = fltarr(nlons-4,nlats-4,naltsprint)
NewnCO = fltarr(nlons-4,nlats-4,naltsprint)
NewnN2 = fltarr(nlons-4,nlats-4,naltsprint)
NewnO2 = fltarr(nlons-4,nlats-4,naltsprint)
NewnO  = fltarr(nlons-4,nlats-4,naltsprint)
;
NewnOP    = fltarr(nlons-4,nlats-4,naltsprint)
NewnO2P   = fltarr(nlons-4,nlats-4,naltsprint)
NewnCO2P  = fltarr(nlons-4,nlats-4,naltsprint)
Newn_e    = fltarr(nlons-4,nlats-4,naltsprint)
;
NewVeast    = fltarr(nlons-4,nlats-4,naltsprint)
NewVnorth   = fltarr(nlons-4,nlats-4,naltsprint)
NewVup      = fltarr(nlons-4,nlats-4,naltsprint)
;
NewQOP1  = fltarr(nlons-4,nlats-4,naltsprint)
NewQO2   = fltarr(nlons-4,nlats-4,naltsprint)
NewQCO2  = fltarr(nlons-4,nlats-4,naltsprint)
NewQN2   = fltarr(nlons-4,nlats-4,naltsprint)
NewQOP2  = fltarr(nlons-4,nlats-4,naltsprint)

;;; We now have a grid in terms of aerodetic
;;; latitude and altitude.
;;; We must interpolate this aerodetic grid
;;; to our desired grid defined by the original
;;; aerocentric grid

for iLon = 0,nlons-5 do begin 
    for iAlt = 0,naltsprint-1 do begin 

        ; set T = lat[*]
        X = reform(AeroGITMLats[*,iAlt])
        T = lat[*]*(180.0/!pi)          

        ;Interpolate our Aerodetic 
        ;latitude grid to
        ;a new more regular grid
        Y = AeroGITMLats[*,iAlt]
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewLats = NewY

        ;Interpolate our Aerodetic 
        ;altitude grid to
        ;a new more regular grid
        Y = AeroGITMAlts[*,iAlt]
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewAlts[*,iAlt] = NewY[*]

        ;;; Non-Logarithmic Splines
        ; interpolate Tn
        Y = reform(Tn[iLon,*,iAlt])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewTn[iLon,0:nlats-5,iAlt] = NewY 

        ; interpolate Te
        Y = reform(Te[iLon,*,iAlt])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewTe[iLon,0:nlats-5,iAlt] = NewY 

        ; interpolate Ti
        Y = reform(Ti[iLon,*,iAlt])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewTi[iLon,0:nlats-5,iAlt] = NewY 

        ; interpolate Veast
        Y = reform(vEast[iLon,*,iAlt])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewVeast[iLon,0:nlats-5,iAlt] = NewY 

        ; interpolate Vnorth
        Y = reform(vNorth[iLon,*,iAlt])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewVnorth[iLon,0:nlats-5,iAlt] = NewY 

        ; interpolate Vnorth
        Y = reform(vUp[iLon,*,iAlt])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewVup[iLon,0:nlats-5,iAlt] = NewY 

        ;;; LOGARITHMIC INTERPOLATION
        ; interpolate nCO2
        Y = reform(alog(nCO2[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewnCO2[iLon,0:nlats-5,iAlt] = exp(NewY) 

        ; interpolate nCO
        Y = reform(alog(nCO[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewnCO[iLon,0:nlats-5,iAlt] = exp(NewY) 
 
        ; interpolate nN2
        Y = reform(alog(nN2[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewnN2[iLon,0:nlats-5,iAlt] = exp(NewY) 

        ; interpolate nO2
        Y = reform(alog(nO2[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewnO2[iLon,0:nlats-5,iAlt] = exp(NewY) 

        ; interpolate nO
        Y = reform(alog(nO[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewnO[iLon,0:nlats-5,iAlt] = exp(NewY) 

        ;;;; Fix for the low ion densities on night side
        SubDen = nOP[iLon,*,iAlt]
        LowDenMask = where(SubDen lt 1.0e-03)
        SubDen[LowDenMask] = 1.0e-03
        nOP[iLon,*,iAlt] = SubDen
        ; interpolate nOP
        Y = reform(alog(nOP[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewnOP[iLon,0:nlats-5,iAlt] = exp(NewY) 

        ;;;; Fix for the low ion densities on night side
        SubDen = nO2P[iLon,*,iAlt]
        LowDenMask = where(SubDen lt 1.0e-03)
        SubDen[LowDenMask] = 1.0e-03
        nO2P[iLon,*,iAlt] = SubDen
        ; interpolate nO2P
        Y = reform(alog(nO2P[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewnO2P[iLon,0:nlats-5,iAlt] = exp(NewY) 

        ;;;; Fix for the low ion densities on night side
        SubDen = nCO2P[iLon,*,iAlt]
        LowDenMask = where(SubDen lt 1.0e-03)
        SubDen[LowDenMask] = 1.0e-03
        nCO2P[iLon,*,iAlt] = SubDen
        ; interpolate nCO2P
        Y = reform(alog(nCO2P[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewnCO2P[iLon,0:nlats-5,iAlt] = exp(NewY) 

        ;;;; Fix for the low ion densities on night side
        SubDen = n_e[iLon,*,iAlt]
        LowDenMask = where(SubDen lt 1.0e-03)
        SubDen[LowDenMask] = 1.0e-03
        n_e[iLon,*,iAlt] = SubDen
        ; interpolate n_e
        Y = reform(alog(n_e[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        Newn_e[iLon,0:nlats-5,iAlt] = exp(NewY) 

;;;; Next, fill in the Ionization Rates if they are missing for
;;;; the nightside.
        SubDen = QOP1[iLon,*,iAlt]
        LowDenMask = where(SubDen lt 1.0e-30)
        SubDen[LowDenMask] = 1.0e-30
        QOP1[iLon,*,iAlt] = SubDen
        ; interpolate QOP1
        Y = reform(alog(QOP1[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewQOP1[iLon,0:nlats-5,iAlt] = exp(NewY) 

        SubDen = QO2[iLon,*,iAlt]
        LowDenMask = where(SubDen lt 1.0e-30)
        SubDen[LowDenMask] = 1.0e-30
        QO2[iLon,*,iAlt] = SubDen
        ; interpolate QO2
        Y = reform(alog(QO2[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewQO2[iLon,0:nlats-5,iAlt] = exp(NewY) 

        SubDen = QCO2[iLon,*,iAlt]
        LowDenMask = where(SubDen lt 1.0e-30)
        SubDen[LowDenMask] = 1.0e-30
        QCO2[iLon,*,iAlt] = SubDen
        ; interpolate QCO2
        Y = reform(alog(QCO2[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewQCO2[iLon,0:nlats-5,iAlt] = exp(NewY) 

        SubDen = QN2[iLon,*,iAlt]
        LowDenMask = where(SubDen lt 1.0e-30)
        SubDen[LowDenMask] = 1.0e-30
        QN2[iLon,*,iAlt] = SubDen
        ; interpolate QN2
        Y = reform(alog(QN2[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewQN2[iLon,0:nlats-5,iAlt] = exp(NewY) 

        SubDen = QOP2[iLon,*,iAlt]
        LowDenMask = where(SubDen lt 1.0e-30)
        SubDen[LowDenMask] = 1.0e-30
        QOP2[iLon,*,iAlt] = SubDen
        ; interpolate QOP2
        Y = reform(alog(QOP2[iLon,*,iAlt]))
        NewY = SPLINE(X,Y,T,/DOUBLE)
        NewQOP2[iLon,0:nlats-5,iAlt] = exp(NewY) 

    endfor 
endfor ;iLon = 0,nlons-5 do begin 


;;;----- Now, all fields have been interpolated to a new
;;;----- Aerodetic Latitude Grid with the same
;;;----- cell values as our old aerocentric latitude grid
     

;-NEXT, We must now take the user's input in altitude
; and generate the desired grid.
ialtmin_user = min(where(alt ge user_alt_min))
ialtmax_user = max(where(alt le user_alt_max))
naltsprint_user = ialtmax_user - ialtmin_user + 1

useraltgrid = alt[ialtmin_user:ialtmax_user]

;;; next define the size of this desired grid:

AerodeticAlts = fltarr(naltsprint_user)
AerodeticLats = fltarr(nlats-4)

AerodeticTn   = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticTe   = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticTi   = fltarr(nlons-4,nlats-4,naltsprint_user)
;
AerodeticVnorth   = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticVeast    = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticVup      = fltarr(nlons-4,nlats-4,naltsprint_user)
;
AerodeticnCO2 = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticnCO  = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticnN2  = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticnO2  = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticnO   = fltarr(nlons-4,nlats-4,naltsprint_user)
;
AerodeticnCO2P = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticnO2P  = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticnOP   = fltarr(nlons-4,nlats-4,naltsprint_user)
Aerodeticn_e   = fltarr(nlons-4,nlats-4,naltsprint_user)
;
AerodeticQOP1  = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticQO2   = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticQCO2  = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticQOP2  = fltarr(nlons-4,nlats-4,naltsprint_user)
AerodeticQN2   = fltarr(nlons-4,nlats-4,naltsprint_user)

;;; Aerodetic Lats is easy, since we already
;;; defined this earlier
AerodeticLats = NewLats
AerodeticAlts = useraltgrid

for iLon = 0,nlons-5 do begin 
    for iLat = 0,nlats-5 do begin 

        X = reform(NewAlts[iLat,*])
        T = useraltgrid

        ;;; NON-LOG-----
        ; Tn
        Y = reform(  NewTn[iLon,iLat,*])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        AerodeticTn[iLon,iLat,0:naltsprint_user-1] = NewY 

        ;Te
        Y = reform(  NewTe[iLon,iLat,*])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        AerodeticTe[iLon,iLat,0:naltsprint_user-1] = NewY 

        ;Ti
        Y = reform(  NewTi[iLon,iLat,*])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        AerodeticTi[iLon,iLat,0:naltsprint_user-1] = NewY 

        ;Veast
        Y = reform(  NewVeast[iLon,iLat,*])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        AerodeticVeast[iLon,iLat,0:naltsprint_user-1] = NewY 

        ;Vnorth
        Y = reform(  NewVnorth[iLon,iLat,*])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        AerodeticVnorth[iLon,iLat,0:naltsprint_user-1] = NewY 

        ;Vup
        Y = reform(  NewVup[iLon,iLat,*])
        NewY = SPLINE(X,Y,T,/DOUBLE)
        AerodeticVup[iLon,iLat,0:naltsprint_user-1] = NewY 


        ;LOGARITHMIC 

        ; nCO2
        Y = reform(  alog(NewnCO2[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticnCO2[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; nCO
        Y = reform(  alog(NewnCO[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticnCO[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; nN2
        Y = reform(  alog(NewnN2[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticnN2[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; nO2
        Y = reform(  alog(NewnO2[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticnO2[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; nO
        Y = reform(  alog(NewnO[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticnO[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; nCO2P
        Y = reform(  alog(NewnCO2P[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticnCO2P[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; nO2P
        Y = reform(  alog(NewnO2P[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticnO2P[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; nOP
        Y = reform(  alog(NewnOP[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticnOP[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; n_e
        Y = reform(  alog(Newn_e[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        Aerodeticn_e[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; QO
        Y = reform(  alog(NewQOP1[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticQOP1[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; QO2
        Y = reform(  alog(NewQO2[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticQO2[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; QCO2
        Y = reform(  alog(NewQCO2[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticQCO2[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; QOP2
        Y = reform(  alog(NewQOP2[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticQOP2[iLon,iLat,0:naltsprint_user-1] = NewY 

        ; QN2
        Y = reform(  alog(NewQN2[iLon,iLat,*]))
        NewY = exp(SPLINE(X,Y,T,/DOUBLE))
        AerodeticQN2[iLon,iLat,0:naltsprint_user-1] = NewY 

    endfor 
endfor ;iLon = 0,nlons-5 do begin 


;mindata = 100.0
;maxdata = 195.0
;
;plotvar = Tn[16,*,*]
;
;levels = 45
;step = (maxdata - mindata)/levels
;UserLevels = indgen(levels+1)*step + mindata
;numlevels = n_elements(UserLevels)
;cgloadct, 74, ncolors = levels + 4, bottom = 2, /SILENT, /REVERSE
;
;clevels = intarr(numlevels+2)
;clevels[1:numlevels] = UserLevels[0:numlevels-1]
;clevels[0] = UserLevels[0]
;clevels[numlevels+1] = UserLevels[numlevels-1]
;
;psfile = 'Aerocentric_Temperature.eps'
;cgps_open, filename = psfile, /quiet
;   cgContour, reform(Tn[16,*,ialtmin_user:ialtmax_user]), lat*(180.0/!pi), $
;          alt[ialtmin_user:ialtmax_user], $
;          ystyle = 1, c_colors = indgen(levels+1)+3,$
;          color = 'black', background = 'white', $
;          /FILL, yrange = [user_alt_min, user_alt_max],$
;          title = textoidl('Aerocentric T_{n} (K)'),$
;          ytitle = 'Aerocentric Altitude (km)',$
;          xtitle = 'Aerocentric Latitude (deg)',$
;          Levels = UserLevels
;
;   cgColorBar, ncolors = levels, bottom = 3, divisions = 10,$
;          range = [mindata,maxdata], format = '(f5.1)',$
;          /VERTICAL, /RIGHT, $
;          position = [0.95, 0.125, 1.00, 0.90]
;
;cgps_close
;     
;
;
;plotvar = AerodeticTn[16,*,*]
;
;levels = 45
;step = (maxdata - mindata)/levels
;UserLevels = indgen(levels+1)*step + mindata
;numlevels = n_elements(UserLevels)
;cgloadct, 74, ncolors = levels + 4, bottom = 2, /SILENT, /REVERSE
;
;clevels = intarr(numlevels+2)
;clevels[1:numlevels] = UserLevels[0:numlevels-1]
;clevels[0] = UserLevels[0]
;clevels[numlevels+1] = UserLevels[numlevels-1]
;
;psfile = 'Aerodetic_Temperature.eps'
;cgps_open, filename = psfile, /quiet
;   cgContour, reform(AerodeticTn[16,*,*]), $
;          AerodeticLats*(180.0/!pi),AerodeticAlts, $
;          ystyle = 1, c_colors = indgen(levels+1)+3,$
;          color = 'black', background = 'white', $
;          /FILL, yrange = [user_alt_min, user_alt_max],$
;          title = textoidl('Aerodetic T_{n} (K)'),$
;          ytitle = 'Aerodetic Altitude (km)',$
;          xtitle = 'Aerodetic Latitude (deg)',$
;          Levels = UserLevels
;
;   cgColorBar, ncolors = levels, bottom = 3, divisions = 10,$
;          range = [mindata,maxdata], format = '(f5.1)',$
;          /VERTICAL, /RIGHT, $
;          position = [0.95, 0.125, 1.00, 0.90]
;
;cgps_close


openw,1,'C:\Users\rshar\IDLWorkspace83\Default\Current Data\April 2017\2D Plot Test Area\gitm_'+date+'_userdetic.dat'

printf, 1, '#MGITM Results on '+tostr(year)+'-'+chopr('0'+tostr(month),2)+'-'+chopr('0'+tostr(day),2)+' at '+ $
	chopr('0'+tostr(hour),2)+':'+chopr('0'+tostr(min),2)+':'+chopr('0'+tostr(sec),2) + ' UT.'  
printf, 1, '#Each column contains the following variables at the given longitude, latitude, and altitude.'
printf,1,'#Number of longitude points: '+tostr(nlons-4)
printf,1,'#Number of latitude points: '+tostr(nlats-4)
printf,1,'#Number of altitude points:  '+tostr(nAltsprint_user)
printf,1,'#Units are SI- Densities: #/m3, temperatures: K, wind velocities : m/s. '
printf,1,'#1.Longitude 2.Latitude 3.Altitude 4.Tn 5.Ti 6.Te 7.nCO2 8.nO 9.nN2 10.nCO, 11.nO2, 12.nO2P 13.nOP 14.nCO2P 15.Ne 16.UN 17.VN 18.WN 19.SZA ' 
printf,1, '#Start'

for ialt = 0, naltsprint_user - 1 do begin
   for ilat = 0, nlats - 5 do begin
      for ilon = 0,nlons - 5 do begin


;		printf, 1, lon(ilon)*180/!pi,lat(ilat)*180/!pi,alt(ialtmin+ialt),$
;                        tn(ilon,ilat,ialt),ti(ilon,ilat,ialt),te(ilon,ilat,ialt),$ 
;			nCO2(ilon,ilat,iAlt),$
;                        nO(ilon,ilat,ialt),$
;                        nN2(ilon,ilat,ialt),$
;                        nCO(ilon,ilat,ialt),$
;                        nO2(ilon,ilat,ialt),$
;                        nO2P(iLon,ilat,iAlt),$
;                        nOP(iLon,ilat,iAlt),$
;                        nCO2P(iLon,ilat,iAlt),$
;			n_e(iLon,ilat,iAlt),$
;                        vEast(iLon,ilat,iAlt), vNorth(iLon,ilat,iAlt), vUp(iLon,ilat,iAlt),$
;			format='(18G12.5)'

		printf, 1, lon(ilon)*180/!pi,AerodeticLats(ilat),$
                        AerodeticAlts(ialtmin+ialt),$
                        AerodeticTn(ilon,ilat,ialt),$
                        AerodeticTi(ilon,ilat,ialt),$
                        AerodeticTe(ilon,ilat,ialt),$
			AerodeticnCO2(ilon,ilat,iAlt),$
			AerodeticnO(ilon,ilat,iAlt),$
			AerodeticnN2(ilon,ilat,iAlt),$
			AerodeticnCO(ilon,ilat,iAlt),$
			AerodeticnO2(ilon,ilat,iAlt),$
			AerodeticnO2P(ilon,ilat,iAlt),$
			AerodeticnOP(ilon,ilat,iAlt),$
			AerodeticnCO2P(ilon,ilat,iAlt),$
			Aerodeticn_e(ilon,ilat,iAlt),$
			AerodeticVeast(ilon,ilat,iAlt),$
			AerodeticVnorth(ilon,ilat,iAlt),$
			AerodeticVup(ilon,ilat,iAlt),$
 			SZA(iLon,ilat,0),$
			format='(19G12.5)'

             endfor
     endfor
endfor
close,1




;meta = {year:year,month:month,day:day,$
;        hour:hour,min:min,sec:sec,nAlts:naltsprint,nLons:nLons-4,nlats:nlats-4,$
;       latitude : reform(data(1,0,2:nlats-3,0))*180/!pi, Altitude:reform(data(2,0,0,ialtmin:ialtmax)/1000.0),$
;        Longitude:reform(data(0,2:nlons-3,0,0))*180/!pi,SZA:SZA,LS:LS,LONGSUBSOL:subsolarlon,$
;        coordinate_system:'GEO',altitude_from:'surface',Mars_Radius:3388.25}

meta = {year:year,month:month,day:day,$
        hour:hour,min:min,sec:sec,nAlts:naltsprint_user,nLons:nLons-4,nlats:nlats-4,$
       latitude : lat*180/!pi, Altitude:alt,$
        Longitude:lon*180/!pi,SZA:SZA,LS:LS,LONGSUBSOL:subsolarlon,$
        coordinate_system:'Aerodetic',altitude_from:'surface',Mars_Radius:3388.25}
;
;NDensityS = {CO2:nCO2,CO:nCO,N2:nN2,O2:nO2,O:nO}
;IDensityS = {O2P:nO2P,OP:nOP,CO2P:nCO2P,n_e:n_e}
;NVelocity =  {VEast:VEast,VNorth:Vnorth,VUp:Vup}
;Temperature = {Tn:Tn,Ti:Ti,Te:Te}
;QEUVIonRate =  {O:QOP1,O2:QOP2,N2:QN2,CO2:QCO2}

NDensityS = {CO2:AerodeticnCO2,CO:AerodeticnCO,N2:AerodeticnN2,O2:AerodeticnO2,O:AerodeticnO}
IDensityS = {O2P:AerodeticnO2P,OP:AerodeticnOP,CO2P:AerodeticnCO2P,n_e:Aerodeticn_e}
NVelocity =  {VEast:AerodeticVEast,VNorth:AerodeticVnorth,VUp:AerodeticVup}
Temperature = {Tn:AerodeticTn,Ti:AerodeticTi,Te:AerodeticTe}
QEUVIonRate =  {O:AerodeticQOP1,O2:AerodeticQOP2,N2:AerodeticQN2,CO2:AerodeticQCO2}


save_filename = 'C:\Users\rshar\IDLWorkspace83\Default\Current Data\April 2017\2D Plot Test Area\gitm_'+date+'_userdetic.sav'
save,meta,Ndensitys,iDensityS,NVelocity,Temperature,QEUVIonRate,filename=save_filename

end
