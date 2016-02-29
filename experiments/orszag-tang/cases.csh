@ c = 10
while ($c < 18)
  ./dispatch.x case$c.nml > case$c.log &
  @ c++
  ./dispatch.x case$c.nml > case$c.log &
  @ c++
  wait
end
