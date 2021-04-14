
if(F)
{
  # wdk
  ml = 500
  vol = .4

  # wdk
  ml = 250
  vol = .4

  # wdk
  ml = 3/4*250
  vol = .4
  
  # pwo
  ml = 500
  vol = .04
  
  # wdk
  ml = 50
  vol = .4
  
  # 4x pwo 0.04
  ml = 500 * 4
  vol = .04
  
  
  # 4x pwo 0.05
  ml = 500 * 4
  vol = .05
}

ml = 50 * 4
vol = .40

start = 19

{
  W = 82
  g = .8
  K = .7
  
  P = ml * vol * g / K / W
  P = round(P, 2)
  h = P / .12
  h = round(h, 1)
  
  cat(sprintf("P = %g \n", P))
  cat(sprintf("h = %g \n", h))
  cat(sprintf("end = %g \n", (start + h) %% 24))
  
}



