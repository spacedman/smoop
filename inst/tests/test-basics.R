
test_that(
  "basic things",
  {
    set.seed(310366)
    d = data.frame(x=runif(10),y=runif(10),
      Count=as.integer(runif(10)*20)
      )
    d$N = 2*max(d$Count)
    coordinates(d)=~x+y
    expect_error(smoop())
    expect_error(smoop(~Count,~N))
    expect_error(smoop(~Count,~N,d)) # M missing
    smoop(~Count,~N,d,max(d$N)*10,nx=10,ny=10)
  }
  )

