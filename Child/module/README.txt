
# for your .cshrc

setenv EXPLORERUSERHOME ${HOME}/exploreruserhome
setenv EXPLORERHOME /usr/local/explorer
if ($?LM_LICENSE_FILE) then
  setenv LM_LICENSE_FILE ${LM_LICENSE_FILE}:/opt/explorer/license/license.dat
else
  setenv LM_LICENSE_FILE /opt/explorer/license/license.dat
endif

if ( $?EXPLORERHOME ) then
  setenv PATH ${EXPLORERHOME}/bin:${PATH}
  # if your X server doesn't support GLX.
  setenv CXGLTYPE mesagl 
endif
