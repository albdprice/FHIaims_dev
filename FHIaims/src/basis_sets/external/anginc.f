C* Module Name: anginc
c
c Numerical Gauss integrations over the unit sphere by a rule
c with octahedral symmetry.
c
c Author: Bernard Delley, Paul Scherrer Institut, Switzerland 
c title = {High Order Integration Schemes on the Unit Sphere}
c Journal of Computational Chemistry 1996, 17, 1152--1155
c
c VB 10/14/05: This subroutine was kindly provided by Bernard Delley
c under the condition that we cite the above work (a matter
c of course, naturally), but no other conditions imposed. 
c
c Input variables:
c  lmax   order of the numerical integration rule
c  ndim   maximum number of integration points
c         (ndim=1202 is sufficient for all schemes here)
c
c Output variables:
c  rpt    integration points
c  wpt    integration weights
c  npt    number of integration points
c           returns -1 if no sufficiently high order available
c                   -n if npt > ndim, n = required ndim
c
      subroutine anginc(rpt,   wpt,   npt,   ndim,    lmax)
 

      implicit double precision (a-h,o-z)
      
      integer npt, ndim, lmax

      double precision    rpt(3,ndim),       wpt(ndim)

      integer   nk5(8),   nk3(8),
     |          nk7(8),
     |          nk11(8),
     1          nk17(8),
     2          nk23(8),
     3          nk29(8),
     4          nk35(8),
     5          nk41(8),
     6          nk47(8),
     7          nk53(8),
     8          nk59(8)

      double precision    pa5(2),  pa3(1),
     |                    pa7(3),
     |                    pa11(5),
     1                    pa17(10),
     2                    pa23(16),
     3                    pa35(33),
     4                    pa29(24),
     5                    pa41(44),
     6                    pa47(56),
     7                    pa53(70),
     8                    pa59(85)

      data nk3,pa3 /  3,   1,   0,   0,   0,   0,   0,    6,
     I  1.6666666666666667D-01 /
      data nk5,pa5 /  5,   1,   1,   0,   0,   0,   0,   14,
     I  6.6666666666666667D-02,  7.5000000000000000D-02 /
      data nk7,pa7 /  7,   1,   1,   1,   0,   0,   0,   26,
     I  4.7619047619047619D-02,  3.2142857142857143D-02,
     I  3.8095238095238095D-02  /
      data nk11,pa11 / 11,   1,   1,   1,   1,   0,   0,  50,
     I  1.2698412698412698D-02,  2.1093750000000000D-02,
     I  2.2574955908289242D-02,  3.0151134457776362D-01,
     I  2.0173335537918871D-02 /
      data nk17,pa17 /   17,   1,   1,   0,   3,   1,   0,  110,
     I  3.8282704949371616D-03,  9.7937375124875125D-03,
     I  1.8511563534473617D-01,  8.2117372831911110D-03,
     I  3.9568947305594191D-01,  9.5954713360709628D-03,
     I  6.9042104838229218D-01,  9.9428148911781033D-03,
     I  4.7836902881215020D-01,  9.6949963616630283D-03  /
      data nk23,pa23 /   23,   1,   1,   1,   4,   1,   1,  194,
     I  1.7823404472446112D-03,  5.5733831788487380D-03,
     I  5.7169059499771019D-03,  4.4469331787174373D-01,
     I  5.5187714672736137D-03,  2.8924656275754386D-01,
     I  5.1582377118053831D-03,  6.7129734426952263D-01,
     I  5.6087040825879968D-03,  1.2993354476500669D-01,
     I  4.1067770281693941D-03,  3.4577021976112827D-01,
     I  5.0518460646148085D-03,  1.5904171053835295D-01,
     I  5.2511857244364202D-01,  5.5302489162330937D-03  /
      data nk29,pa29 /   29,   1,   1,   0,   6,   2,   2,  302,
     I  8.5459117251281481D-04,  3.5991192850255715D-03,
     I  7.0117664160895449D-01,  3.6500458076772554D-03,
     I  6.5663294102196118D-01,  3.6048226014198817D-03,
     I  4.7290541325810046D-01,  3.5767296617433671D-03,
     I  3.5156403455701051D-01,  3.4497884243058833D-03,
     I  2.2196452362941784D-01,  3.1089531224136753D-03,
     I  9.6183085226147838D-02,  2.3521014136891644D-03,
     I  5.7189558918789607D-01,  3.6008209322164603D-03,
     I  2.6441528870606625D-01,  2.9823449631718039D-03,
     I  2.5100347517704651D-01,  5.4486773725807738D-01,
     I  3.5715405542733871D-03,  1.2335485325833274D-01,
     I  4.1277240831685310D-01,  3.3923122050061702D-03  /
      data nk35,pa35 /   35,   1,   1,   1,   7,   2,   4,  434,
     I  5.2658979682244362D-04,  2.5123174189273072D-03,
     I  2.5482199720026072D-03,  6.9093463075091106D-01,
     I  2.5304038011863550D-03,  6.4566647074242561D-01,
     I  2.5132671745975644D-03,  4.9143426377847465D-01,
     I  2.5017251684029361D-03,  3.9272597633680022D-01,
     I  2.4453734373129800D-03,  2.8612890103076384D-01,
     I  2.3026947822274158D-03,  1.7748360546091578D-01,
     I  2.0142790209185282D-03,  7.5680843671780184D-02,
     I  1.4624956215946138D-03,  2.1027252285730696D-01,
     I  1.9109512821795323D-03,  4.7159869115131592D-01,
     I  2.4174423756389808D-03,  9.9217696364292373D-02,
     I  3.3443631453434549D-01,  2.2366077604378487D-03,
     I  2.0548236964030437D-01,  4.5023303825826254D-01,
     I  2.4169300443247753D-03,  3.1042840351665415D-01,
     I  5.5501523610768072D-01,  2.4966440545530860D-03,
     I  1.0680182607580483D-01,  5.9051570489252711D-01,
     I  2.5122368545634951D-03  /
c alternate solution 950704 BD
      data nk41,pa41 /   41,   1,   1,   0,   9,   3,   6,   590,
     &  3.0951212953061873D-04,  1.8523796985974890D-03,
     &  6.0950341155071959D-02,  9.7643311650510500D-04,
     &  1.4590364491577632D-01,  1.3847372348516919D-03,
     &  2.3847367014218874D-01,  1.6172106472544112D-03,
     &  3.3179207364721231D-01,  1.7495646572811541D-03,
     &  4.2157617840109665D-01,  1.8184717781627688D-03,
     &  5.0444197078003583D-01,  1.8467159561512418D-03,
     &  6.3725469392587524D-01,  1.8520288282962131D-03,
     &  6.8077440664552429D-01,  1.8588125854383170D-03,
     &  7.0409549382274691D-01,  1.8717906392777438D-03,
     &  1.7247820099077235D-01,  1.3003216858860477D-03,
     &  3.9647553481998576D-01,  1.7051539963958640D-03,
     &  6.1168434420098755D-01,  1.8571611967740780D-03,
     &  8.2130215819325114D-02,  2.7786731905862443D-01,
     &  1.5552136033968085D-03,  8.9992058420748749D-02,
     &  5.0335642710751172D-01,  1.8022391280085255D-03,
     &  1.8166408403602095D-01,  5.9841264978853796D-01,
     &  1.8498305604436602D-03,  1.7207952256568781D-01,
     &  3.7910354076955633D-01,  1.7139045071067087D-03,
     &  2.6347166559379496D-01,  4.7423928425519802D-01,
     &  1.8026589343774512D-03,  3.5182809277335190D-01,
     &  5.6102638086220602D-01,  1.8428664729052856D-03  /
c alternate solution 950705 BD
      data nk47,pa47 /   47,   1,   1,   1,  10,   3,   9,   770,
     &  2.1929420881811841D-04,  1.4219403443358774D-03,
     &  1.4364336173190798D-03,  5.0872044105023605D-02,
     &  6.7981235110505020D-04,  1.2281987901788307D-01,
     &  9.9131842352949122D-04,  2.0268908144087861D-01,
     &  1.1802078332389488D-03,  2.8477451564642939D-01,
     &  1.2965996020809207D-03,  3.6567190789780265D-01,
     &  1.3658714274283164D-03,  4.4282648867134686D-01,
     &  1.4029886047753253D-03,  5.1406196272497354D-01,
     &  1.4186455635956094D-03,  6.3064012191668026D-01,
     &  1.4213767418516618D-03,  6.7168833320226119D-01,
     &  1.4239964754909616D-03,  6.9797926853368807D-01,
     &  1.4315540421785668D-03,  1.4468656741953093D-01,
     &  9.2544014998653679D-04,  3.3902634754112157D-01,
     &  1.2502399950535093D-03,  5.3358046512635063D-01,
     &  1.3943658433292301D-03,  6.9440243933494130D-02,
     &  2.3551878942423264D-01,  1.1270890946717488D-03,
     &  2.2690041095294599D-01,  4.1021824740457302D-01,
     &  1.3457537609106701D-03,  8.0255746077753389D-02,
     &  6.2143024174816046D-01,  1.4249572833167828D-03,
     &  1.4679995278965720D-01,  3.2452843457173944D-01,
     &  1.2615233412377500D-03,  1.5715077698247271D-01,
     &  5.2244821896966297D-01,  1.3925471060526959D-03,
     &  2.3657029931572456D-01,  6.0175466340895581D-01,
     &  1.4187616778776564D-03,  7.7148158667657320D-02,
     &  4.3465755161411628D-01,  1.3383666844795541D-03,
     &  3.0629366662107302D-01,  4.9088265890376162D-01,
     &  1.3937008626761314D-03,  3.8224773795247870D-01,
     &  5.6487681490995005D-01,  1.4159147574669320D-03  /
      data nk53,pa53 /   53,   1,   1,   0,  12,   4,  12,  974,
     I  1.4382941905274311D-04,  1.1257722882870041D-03,
     I  4.2929635453413471D-02,  4.9480293419492410D-04,
     I  1.0514268540864042D-01,  7.3579901091254705D-04,
     I  1.7500248676230874D-01,  8.8891327713043843D-04,
     I  2.4776533796502568D-01,  9.8883478389214349D-04,
     I  3.2065671239559574D-01,  1.0532996817094706D-03,
     I  3.9165207498499835D-01,  1.0927788070145785D-03,
     I  4.5908258741876237D-01,  1.1143893940632272D-03,
     I  5.2145638884158605D-01,  1.1237247880515553D-03,
     I  6.2531702446541989D-01,  1.1252393252438136D-03,
     I  6.6379267445231699D-01,  1.1261532718159050D-03,
     I  6.9104103984983007D-01,  1.1302869311238408D-03,
     I  7.0529070074577603D-01,  1.1349865343639549D-03,
     I  1.2366867626579899D-01,  6.8233679271099310D-04,
     I  2.9407771144683870D-01,  9.4541581604470958D-04,
     I  4.6977538492076491D-01,  1.0744299753856791D-03,
     I  6.3345632411395669D-01,  1.1293000865691317D-03,
     I  5.9740486141813418D-02,  2.0291287527775228D-01,
     I  8.4368845009019544D-04,  1.3757604084736365D-01,
     I  4.6026219424840539D-01,  1.0752557204488846D-03,
     I  3.3910165263362857D-01,  5.0306739996620357D-01,
     I  1.1085772368644620D-03,  1.2716751914398195D-01,
     I  2.8176064224421343D-01,  9.5664753237833573D-04,
     I  2.6931207404135125D-01,  4.3315612917201574D-01,
     I  1.0806632507173907D-03,  1.4197864526019183D-01,
     I  6.2561673585808142D-01,  1.1267971311962946D-03,
     I  6.7092846007382550D-02,  3.7983952168591567D-01,
     I  1.0225687153580612D-03,  7.0577381832561723D-02,
     I  5.5175054214235205D-01,  1.1089602677131075D-03,
     I  2.7838884778821546D-01,  6.0296191561591869D-01,
     I  1.1227906534357658D-03,  1.9795789389174069D-01,
     I  3.5896063295890958D-01,  1.0324018471174598D-03,
     I  2.0873070611032740D-01,  5.3486664381354765D-01,
     I  1.1072493822838539D-03,  4.0551221378728359D-01,
     I  5.6749975460743735D-01,  1.1217800485199721D-03  /
      data nk59,pa59 /   59,   1,   1,   1,  13,   4,  16, 1202,
     I  1.1051892332675715D-04,  9.1331597864435614D-04,
     I  9.2052327380907415D-04,  3.7126364496570891D-02,
     I  3.6904218980178990D-04,  9.1400604122622234D-02,
     I  5.6039909286806603D-04,  1.5310778524699062D-01,
     I  6.8652976292826086D-04,  2.1809288916606116D-01,
     I  7.7203385511456304D-04,  2.8398745322001746D-01,
     I  8.3015459588947951D-04,  3.4911776009637644D-01,
     I  8.6866925501796284D-04,  4.1214314614443092D-01,
     I  8.9270762858468901D-04,  4.7189936271491266D-01,
     I  9.0608202385682188D-04,  5.2731454528423366D-01,
     I  9.1197772549408672D-04,  6.2094753324440192D-01,
     I  9.1287201386041811D-04,  6.5697227118572905D-01,
     I  9.1307149356917351D-04,  6.8417883090701434D-01,
     I  9.1528737845541164D-04,  7.0126043301236308D-01,
     I  9.1874362743216541D-04,  1.0723822154781661D-01,
     I  5.1769773129656942D-04,  2.5820689594969680D-01,
     I  7.3311436821014169D-04,  4.1727529553067168D-01,
     I  8.4632328363799285D-04,  5.7003669117925033D-01,
     I  9.0311226942539918D-04,  5.2106394770112841D-02,
     I  1.7717740226153253D-01,  6.4857784531632566D-04,
     I  1.1156409571564867D-01,  2.4757164634262876D-01,
     I  7.4350309109823692D-04,  1.7465516775786261D-01,
     I  3.1736152466119767D-01,  8.1017314974680177D-04,
     I  2.3902784793817240D-01,  3.8542911506692237D-01,
     I  8.5562992573118124D-04,  3.0294669735289819D-01,
     I  4.5074225931570644D-01,  8.8502823412654443D-04,
     I  3.6498322605976536D-01,  5.1235184864198708D-01,
     I  9.0226929384269151D-04,  4.2386447815223403D-01,
     I  5.6937024984684411D-01,  9.1057602589701256D-04,
     I  5.9058888532355084D-02,  3.3546162890664885D-01,
     I  7.9985278918390537D-04,  1.2172350510959870D-01,
     I  4.0902684270853572D-01,  8.4833895745943309D-04,
     I  1.8575051945473351D-01,  4.7853206759224352D-01,
     I  8.8110481824257202D-04,  2.4941121623622365D-01,
     I  5.4343035696939004D-01,  9.0100916771050857D-04,
     I  3.1122759471496082D-01,  6.0311616930963100D-01,
     I  9.1078135794827047D-04,  6.2662506241541695D-02,
     I  4.9322211848512846D-01,  8.8032086797382601D-04,
     I  1.2677748006842827D-01,  5.6321230207620997D-01,
     I  9.0213422990406534D-04,  1.9060182227792370D-01,
     I  6.2698055090243917D-01,  9.1315780031894351D-04,
     I  6.4245492242205886D-02,  6.3942796347491023D-01,
     I  9.1580161746934653D-04  /

      if(ndim .le. 0) stop 'ndim<=0'
      if(lmax .le. 3) then
         npt = nk3(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk3,  pa3,   rpt,   wpt)
         endif

      else if(lmax .le.  5) then
         npt = nk5(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk5,  pa5,   rpt,   wpt)
         endif

      else if(lmax .le.  7) then
         npt = nk7(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk7,  pa7,   rpt,   wpt)
         endif

      else if(lmax .le. 11) then
         npt = nk11(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk11,  pa11,   rpt,   wpt)
         endif

      else if(lmax .le. 17) then
         npt = nk17(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk17,  pa17,   rpt,   wpt)
         endif

      else if(lmax .le. 23) then
         npt = nk23(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk23,  pa23,   rpt,   wpt)
         endif

      else if(lmax .le. 29) then
         npt = nk29(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk29,  pa29,   rpt,   wpt)
         endif

      else if(lmax .le. 35) then
         npt = nk35(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk35,  pa35,   rpt,   wpt)
         endif

      else if(lmax .le. 41) then
         npt = nk41(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk41,  pa41,   rpt,   wpt)
         endif

      else if(lmax .le. 47) then
         npt = nk47(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk47,  pa47,   rpt,   wpt)
         endif

      else if(lmax .le. 53) then
         npt = nk53(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk53,  pa53,   rpt,   wpt)
         endif

      else if(lmax .le. 59) then
         npt = nk59(8)
         if(npt .gt.ndim) then
            npt = -npt
         else
            call anmesh(   nk59,  pa59,   rpt,   wpt)
         endif

      else
c        such high order is not available here
         npt = -1

      endif

      return
      end
c
C* Module Name: anmesh
c
c  generates integration mesh for integrations over the unit sphere
c  with octahedral symmetry.
c
c Author: Bernard Delley, Paul Scherrer Institut, Switzerland 
c title = {High Order Integration Schemes on the Unit Sphere}
c Journal of Computational Chemistry 1996, 17, 1152--1155

c  Input:
c   nk(8)  vector defining type of integration rule
c   pa(*)  vector defining parameters of integration rule
c
c  Output:
c   rpt(3,*)   cartesian points on the unit sphere
c   wpt(*)     integration weight factors

      subroutine anmesh(  nk,     pa,    rpt,     wpt)

      implicit real*8(a-h,o-z),integer(i-n)

      integer   nk(8) 

      double precision    pa(*),   rpt(3,*),    wpt(*)

      data zero,one,two,three / 0.0d0, 1.0d0, 2.0d0, 3.0d0 /

c     write(6,*) nk
      i = 0
      ip = 0
      if(nk(2) .gt. 0) then   ! 6 points
         ip = ip + 1
         do ix=1,3
            do iy=1,-1,-2
               i = i + 1
               wpt(i) = pa(ip)
               do j=1,3
                  rpt(j,i) = zero
               enddo
               rpt(ix,i) = dfloat(iy)
            enddo
         enddo
      endif

      if(nk(3) .gt. 0) then   ! 8 points
         c = one/dsqrt(three)
         ip = ip + 1
         do ix=1,-1,-2
            do iy=1,-1,-2
               do iz=1,-1,-2
                  i = i + 1
                  wpt(i) = pa(ip)
                  rpt(1,i) = dfloat(ix)*c
                  rpt(2,i) = dfloat(iy)*c
                  rpt(3,i) = dfloat(iz)*c
               enddo
            enddo
         enddo
      endif

      if(nk(4) .gt. 0) then   ! 12 points
         c = one/dsqrt(two)
         ip = ip + 1
         do ix=1,-1,-2
            do iy=1,-1,-2
               do iz=1,3
                  i = i + 1
                  wpt(i) = pa(ip)
                  rpt(iz,i) = dfloat(ix)*c
                  j = mod(iz,3) + 1
                  rpt(j,i) = dfloat(iy)*c
                  j = 6 - iz - j
                  rpt(j,i) = zero
               enddo
            enddo
         enddo
      endif

      n1 = nk(5)              ! 24a points
      do jj=1,n1
         ip = ip + 1
         uu = pa(ip)
         vv = dsqrt(one - two*uu*uu)  ! data were checked
         ip = ip + 1
         do ix=1,-1,-2
            do iy=1,-1,-2
               do iz=1,-1,-2
                  do j=1,3
                     i = i + 1
                     wpt(i) = pa(ip)
                     do j1=1,3
                        rpt(j1,i) = uu
                     enddo
                     rpt(j,i) = vv
                     rpt(1,i) = rpt(1,i)*dfloat(ix)
                     rpt(2,i) = rpt(2,i)*dfloat(iy)
                     rpt(3,i) = rpt(3,i)*dfloat(iz)
                  enddo
               enddo
            enddo
         enddo
      enddo

      n1 = nk(6)              ! 24b points
      do jj=1,n1
         ip = ip + 1
         pp = pa(ip)
         qq = dsqrt(one - pp*pp)  ! data were checked
         ip = ip + 1
         do ix=1,-1,-2
            do iy=1,-1,-2
               do ii=0,1
                  do j=1,3
                     i = i + 1
                     wpt(i) = pa(ip)
                     j1 = mod(j+ii,3)+1
                     rpt(j1,i) = pp*dfloat(ix)
                     j1 = mod(j+1-ii,3) + 1
                     rpt(j1,i) = qq*dfloat(iy)
                     rpt(j,i) = zero
                  enddo
               enddo
            enddo
         enddo
      enddo

      n1 = nk(7)              ! 48 points
      do jj=1,n1
         ip = ip + 1
         rr = pa(ip)
         ip = ip + 1
         ss = pa(ip)
         tt = dsqrt(one - rr*rr - ss*ss)
         ip = ip + 1
         do ix=1,-1,-2
            do iy=1,-1,-2
               do iz=1,-1,-2
                  do j=1,3
                     do ii=0,1
                        i = i + 1
                        wpt(i) = pa(ip)
                        rpt(j,i) = rr*dfloat(ix)
                        j1 = mod(j+ii,3) + 1
                        rpt(j1,i) = ss*dfloat(iy)
                        j1 = mod(j+1-ii,3) + 1
                        rpt(j1,i) = tt*dfloat(iz)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
