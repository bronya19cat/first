c three kinds of way to calculate run-time
cccccccccccccccccccccccccccccccccccccccccccccccccc
c character
c date and time
      call date_and_time(date,time(1))
c body
      call date_and_time(date,time(2))
      read(time(1),*) time1
      read(time(2),*) time2
      print *,date,time(1),"/",time(2),time2-time1
ccccccccccccccccccccccccccccccccccccccccccccccccc
c real
c unit is s
c ? t1==t2
      call cpu_time(t1)
c body
      call cpu_time(t2)
      print '("cpu_time",f6.3,f6.3,f6.3)',t1,t2,t2-t1
cccccccccccccccccccccccccccccccccccccccccccccccccc
c integer
c unit is ms
      call system_clock(it1)
c body
      call system_clock(it2)
      print *,"sysClock",it1,it2,it2-it1