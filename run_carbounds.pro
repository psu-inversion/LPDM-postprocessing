pro run_carbounds,year_tag=year_tag,month_tag=month_tag,run_tag=run_tag

; Total number of towers
n_towers=5

; Total number of groups
n_groups=1

info2={tower_names:strarr(n_towers),lon:fltarr(n_towers),lat:fltarr(n_towers),$
      alt:fltarr(n_towers),outdir:' ',indir:' ',year:0l,month:0l,$
      max_towers:fltarr(n_groups),lag_towers:0l,lpdm_timestep:0d,dx:0d,$
      length:0l,dimx:0d,dimy:0d,dimbx:0d,dimby:0d,dimbz:0d,lag:0d,$
      num_days:0l,wrf_file:' ',rel_rate:0d,num_file_per_h:0d,run:' ',$
      proj_std_lon:0.,proj_res:30e3,hemisphere:1,true_lats:fltarr(2),proj_proj:1,$
      flux_window:0l}
print,tag_names(info2)

restore,"lpd_config.sav",/RELAXED_STRUCTURE_ASSIGNMENT,/VERBOSE,/NO_COMPILE
STRUCT_ASSIGN,/VERBOSE,info,info2
info=info2
print,tag_names(info)

info.run=run_tag
; carsurf input is LPDM output
info.indir=info.outdir
;info.indir='/mc1s2/s3/tul5/OUTPUTS_LPDM_ACT/'
info.outdir='/mc1s2/s4/dfw5129/data/LPDM_2010_fpbounds/'+info.run+'/'
print,info.indir

info.year=year_tag
info.month=month_tag
; Indicate the number of towers in each group
; read in from saved LPDM configuration
;info.max_towers=[1] ;,4,5]  

; LPDM output time step (in seconds)
; int(60/info.num_file_per_h)
; Keep in sync with updates to run_lpdm.pro
; TODO: change to read LPD config file
info.lpdm_timestep=60.

; Resolution of the domain
; km?
info.dx=27.

; Number of particles released over an hour
info.rel_rate=35.

; Number of LPDM output files per hour
; copy from LPDM dtoutp
info.num_file_per_h=3.

; WRF file for projection parameters
info.wrf_file='/mc1s2/s4/mpb186/CMS_WRF/Output/2010012312/wrfout_d01_2010-01-24_12:00:00'

; Lag for each period (spin-up time between flux steps and obs time)
; (in day)
;; flux lag in other scripts
info.lag=info.flux_window

; Number of days
day_month=[31,28,31,30,31,30,31,31,30,31,30,31]
test_leap=leap_or_not(info.year)
if test_leap eq 1 then day_month=[31,29,31,30,31,30,31,31,30,31,30,31]

info.num_days=(day_month(info.month-1)+info.lag)

; Flux time steps (in hour)
; number of hours back to include in the footprint files
;; Make sure this corresponds to the lag
;; carsurf_loop divides this by the flux interval to get dimension
;; sizes.
info.length=info.lag * 24

; Names of the sites (now read in from LPDM config savefile)
;; info.site_names=[strtrim(fix(findgen(n_towers))+1,1)]
;; for i=0,n_towers-1 do if ((i+1) lt 10) then info.site_names[i]='0'+info.site_names[i]

; Copy these from LPDM namelist specification in run_lpdm.pro
; Dimensions of the domain
info.dimx=249
info.dimy=184

; for output of carbound
; Dimensions of the boundaries
; Number of walls
info.dimbx=4.
; Number of horiz. pixels on each boundary
info.dimby=29.
; Number of vertical levels
info.dimbz=6.

info.lag_towers=0l

for i=0,n_groups-1 do begin
    prep_carbound,i,info,input

    create_sub_surf,i,info,input

    submit_carsurf,i,info,input

;    create_sub_bound,i,info,input

;    submit_carbound,i,info,input

    info.lag_towers+=info.max_towers(i)

    print,input.site_names
end


stop
end

pro submit_carbound,group,info,input

spawn,'qsub -q kjd1 Sub_Carb_'+info.run+'_'+strtrim(info.year,1)+'_M'+input.month+'_G'+strtrim(group+1,1)+'.csh'

end

pro create_sub_bound,group,info,input

month=strtrim( fix(  ((findgen(14)+11) mod 12)+1),1)
for i=0,13 do if ((i+11) mod 12)+1 lt 10 then month(i)='0'+month(i)

for j=0,n_elements(info.max_towers)-1 do begin

    openw,lun,'Sub_Carb_'+info.run+'_'+strtrim(info.year,1)+'_M'+month(info.month)+'_G'+strtrim(j+1,1)+'.csh',/get_lun
    printf,lun,"#!/bin/tcsh"
    printf,lun,"#"
    printf,lun,"# Example PBS script to run a job on the mc1 cluster."
    printf,lun,"# The lines beginning #PBS set various queuing parameters."
    printf,lun,"#"
    printf,lun,"#PBS -N CARB_"+info.run+"_"+strtrim(info.year,1)+"_M"+month(info.month)+"_G"+strtrim(j+1,1)
    printf,lun,"#PBS -e /home/mc2/dfw5129/LPDM/logCarbound"+strtrim(info.year,1)+"_M"+month(info.month)+"_G"+strtrim(j+1,1)+".e"
    printf,lun,"#PBS -o /home/mc2/dfw5129/LPDM/logCarbound"+strtrim(info.year,1)+"_M"+month(info.month)+"_G"+strtrim(j+1,1)+".o"
    printf,lun,"#PBS -l nodes=1:ppn=8"
    printf,lun,"#PBS -m bea"
    printf,lun,"#"
    printf,lun,"# o Export all my environment variables to the job"
    printf,lun,"#PBS -V"
    printf,lun,"#"
    printf,lun,"cd /home/mc2/dfw5129/LPDM/"
    printf,lun,"idl << EOF"
    printf,lun,".comp carbounds_loop"
    printf,lun,".comp carbounds_loop"
    printf,lun,"carbounds_loop, tag='"+info.outdir+"input_"+strtrim(info.year,1)+"_M"+input.month+"_G"+strtrim(group+1,1)+".sav'"
    printf,lun,"EOF"
    free_lun,lun
end


end


pro submit_carsurf,group,info,input

spawn,'qsub -q kjd1 Sub_Cars_'+info.run+'_'+strtrim(info.year,1)+'_M'+input.month+'_G'+strtrim(group+1,1)+'.csh'

end

pro create_sub_surf,group,info,input

month=strtrim( fix(  ((findgen(14)+11) mod 12)+1),1)
for i=0,13 do if ((i+11) mod 12)+1 lt 10 then month(i)='0'+month(i)

for j=0,n_elements(info.max_towers)-1 do begin

    openw,lun,'Sub_Cars_'+info.run+'_'+strtrim(info.year,1)+'_M'+month(info.month)+'_G'+strtrim(j+1,1)+'.csh',/get_lun
    printf,lun,"#!/bin/tcsh"
    printf,lun,"#"
    printf,lun,"# Example PBS script to run a job on the mc1 cluster."
    printf,lun,"# The lines beginning #PBS set various queuing parameters."
    printf,lun,"#"
    printf,lun,"#PBS -N CARS_"+info.run+"_"+strtrim(info.year,1)+"_M"+month(info.month)+"_G"+strtrim(j+1,1)
    printf,lun,"#PBS -e /home/mc2/dfw5129/LPDM/logCarsurf"+strtrim(info.year,1)+"_M"+month(info.month)+"_G"+strtrim(j+1,1)+".e"
    printf,lun,"#PBS -o /home/mc2/dfw5129/LPDM/logCarsurf"+strtrim(info.year,1)+"_M"+month(info.month)+"_G"+strtrim(j+1,1)+".o"
    ;; Python:
    ;; 2 week lag with 35 particles/time step, 20s time step, 27 km grid uses about 6GB and 20 hours (probably)
    ;; no idea how to generalize that.
    printf,lun,"#PBS -l nodes=1"
    printf,lun,"#PBS -M dfw5129@psu.edu"
    printf,lun,"#PBS -m bae"
    printf,lun,"#"
    printf,lun,"# o Export all my environment variables to the job"
    printf,lun,"#PBS -V"
    printf,lun,"#"
    printf,lun,"cd /home/mc2/dfw5129/LPDM/"
    printf,lun,". ~/LPDM.intel"
    printf,lun,"module switch python/3.4.0"
    printf,lun,". ~/python34"
    printf,lun,"python3.4 carsurf_loop.py '"+info.outdir+"input_"+strtrim(info.year,1)+"_M"+input.month+"_G"+strtrim(group+1,1)+".sav'"
    free_lun,lun
    print,"Using python3.4"
end


end


pro prep_carbound,group,info,input

here_towers = info.max_towers[group]

input={site_names:strarr(here_towers),lpdm_timestep:0d,$
       dx:0d,length:0l,dimx:0d,dimy:0d,dimbx:0d,dimby:0d,dimbz:0d,$
       lag:0d,num_days:0l,indir:' ',outdir:' ',ini_date:0l,month:' ',$
       wrf_file:' ',rel_rate:0d,num_file_per_h:0d,$
       lon:fltarr(here_towers),lat:fltarr(here_towers),$
       alt:fltarr(here_towers),proj_std_lon:0.,proj_res:30e3,hemisphere:1,$
       true_lats:fltarr(2),proj_proj:1,year:' '}

input.site_names[0:info.max_towers[group]-1]=$
   info.tower_names[info.lag_towers:info.lag_towers+info.max_towers(group)-1]
input.lpdm_timestep=info.lpdm_timestep
input.dx=info.dx
input.length=info.length
input.dimx=info.dimx
input.dimy=info.dimy
input.num_days=info.num_days
input.dimbx=info.dimbx
input.dimby=info.dimby
input.dimbz=info.dimbz
input.rel_rate=info.rel_rate
input.num_file_per_h=info.num_file_per_h
input.lag=info.lag

; copy info needed for aux coords
input.lon[0:info.max_towers[group]-1]=info.lon[info.lag_towers:info.lag_towers+$
                                               info.max_towers(group)-1]
input.lat[0:info.max_towers[group]-1]=info.lat[info.lag_towers:info.lag_towers+$
                                               info.max_towers(group)-1]
input.alt[0:info.max_towers[group]-1]=info.alt[info.lag_towers:info.lag_towers+$
                                               info.max_towers(group)-1]
input.proj_std_lon=info.proj_std_lon
input.proj_res=info.proj_res
input.hemisphere=info.hemisphere
input.true_lats=info.true_lats
input.proj_proj=info.proj_proj
input.year=info.year

spawn,'mkdir -p '+info.outdir+strtrim(info.year,1)
input.wrf_file=info.outdir+strtrim(info.year,1)+'/'+$
               strmid(info.wrf_file,strpos(info.wrf_file,'wrfout_d0'),45)
; I don't need all the file, so I'm dropping everything but XLAT and XLONG
spawn,'ncks -v XLAT,XLONG --dimension time,0 --deflate 6 --netcdf4 '+ $
      info.wrf_file + ' --overwrite --output ' + input.wrf_file

month=strtrim( fix(  ((findgen(14)+11) mod 12)+1),1)
for i=0,13 do if ((i+11) mod 12)+1 lt 10 then month(i)='0'+month(i)
input.month=month(info.month)

input.ini_date=1l
input.outdir=info.outdir+'/'+strtrim(info.year,1)+'/'+input.month+'/GROUP'+strtrim(group+1,1)+'/'
input.indir=info.indir+'/'+strtrim(info.year,1)+'/'+input.month+'/GROUP'+strtrim(group+1,1)+'/'

spawn,'mkdir -p '+info.outdir+strtrim(info.year,1)+'/'+input.month+'/'
spawn,'mkdir -p '+input.outdir+'/'

save,input,filename=info.outdir+'input_'+strtrim(info.year,1)+'_M'+input.month+'_G'+strtrim(group+1,1)+'.sav'

end

function leap_or_not,year_select_in
; This function tests if the year is a leap year or not
; It should work for any year (not limited to a period)
flag_leap=0l

if (((year_select_in mod 4) eq 0) and ((year_select_in mod 100) ne 0)) or $
   ((year_select_in mod 400) eq 0) then flag_leap=1

return, flag_leap
end
