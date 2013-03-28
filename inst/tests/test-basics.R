
test_that(
  "basic things",
  {
    set.seed(310366)
    d = data.frame(x=runif(10),y=runif(10),
      Count=as.integer(runif(10)*20)
      )
    d$XC=d$x
    d$YC=d$y
    d$N = 2*max(d$Count)
    coordinates(d)=~x+y
    expect_error(smoop())
    expect_error(smoop(~Count,~N))
    expect_error(smoop(~Count,~N,d)) # M missing
    s=smoop(~Count,~N,d,max(d$N)*10,nx=10,ny=10)
    expect_equal(sum(values(s[[1]])),983.5537849)
    expect_equal(sum(values(s[[2]])),15.24603904)
    expect_equal(sum(values(s[[3]])),-162.9596449)
### formula manipulations need wrapping in I()
    expect_error(smoop(~XC+YC,~N,d,max(d$N)*10,nx=10,ny=10))
    expect_error(smoop(~XC,~N+Count,d,max(d$N)*10,nx=10,ny=10))
    smoop(~I(XC+YC),~N,d,max(d$N)*10,nx=10,ny=10)
### M toooo big
    expect_error(smoop(~Count,~N,d,sum(d$N)+100,nx=10,ny=10))
  }
  )

