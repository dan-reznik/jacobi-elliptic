!===============================================================================
!  This file contains Fortran 90 programs to compute five complete elliptic
!  integrals of first and second kind: K(m), E(m), B(m), D(m), S(m)
!
!     1) xcei.f90  test driver of "ceik", "ceie", "ceib", "ceid", and "ceis"
!     2) ceik.f90  real*8 function to compute K(m)
!     3) ceie.f90  real*8 function to compute E(m)
!     4) ceib.f90  real*8 function to compute B(m) = (E(m)-(1-m)K(m))/m
!     5) ceid.f90  real*8 function to compute D(m) = (K(m)-E(m))/m
!     6) ceis.f90  real*8 function to compute S(m) = (D(m)-B(m))/m
!
!  Notice that B(m) and S(m) are most fundamental since the others are
!  computable from them without suffering precision loss when m is small
!
!  The file also includes the output file of "xcei" 
!
!===============================================================================
program xcei
!
!   Test driver of complete elliptic integrals of first and second kind:
!   K(m), E(m), B(m), D(m), S(m)
!
!     References:
!
!       T. Fukushima, (2015), J. Comp. Appl. Math., 282, 71-76
!
!      "Precise and Fast Computation of General Complete Elliptic
!       Integral of Second Kind"
!
!       T. Fukushima, (2016), Astron. J., re-revised
!
!      "Zonal Toroidal Harmonic Expansions of External Gravitational
!       Fields for Ring-Like Objects"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
nend=20
dmc=1.d0/dble(nend)
!
write (*,"(a10,5a25)") "m","K(m)","E(m)","B(m)","D(m)","S(m)"
do n=1,nend
    fmc=dmc*dble(n)
    cek=ceik(fmc)
    cee=ceie(fmc)
    ceb=ceib(fmc)
    ced=ceid(fmc)
    ces=ceis(fmc)
    fm=1.d0-fmc
    write (*,"(0pf10.5,1p5e25.15)") fm,cek,cee,ceb,ced,ces
enddo
!
stop
end program xcei
!===============================================================================
real*8 function ceik(mc)
!
! Double precision minimax rational approximation of K(m):
! the complete elliptic integral of the first kind
!
!     Reference: T. Fukushima, (2015), J. Comp. Appl. Math., 282, 71-76
!
!      "Precise and Fast Computation of General Complete Elliptic
!       Integral of Second Kind"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
real*8 mc
real*8 t,t2
if(mc.gt.1.d0) then
    write(*,"(a20,1pe15.7)") "(ceik) negative m: mc=",mc
elseif(mc.gt.0.592990d0) then
	t=2.45694208987494165d0*mc-1.45694208987494165d0
	t2=t*t
    ceik=((3703.75266375099019d0 &
    +t2*(2744.82029097576810d0 &
    +t2*36.2381612593459565d0)) &
    +t*(5462.47093231923466d0 &
    +t2*(543.839017382099411d0 &
    +t2*0.393188651542789784d0)))/ &
    ((2077.94377067058435d0 &
    +t2*(1959.05960044399275d0 &
    +t2*43.5464368440078942d0)) &
    +t*(3398.00069767755460d0 &
    +t2*(472.794455487539279d0 &
    +t2)))
elseif(mc.gt.0.350756d0) then
    t=4.12823963605439369d0*mc-1.44800482178389491d0
    t2=t*t
    ceik=((4264.28203103974630d0 &
    +t2*(3214.59187442783167d0 &
    +t2*43.2589626155454993d0)) &
    +t*(6341.90978213264024d0 &
    +t2*(642.790566685354573d0 &
    +t2*0.475223892294445943d0)))/ &
    ((2125.06914237062279d0 &
    +t2*(2006.03187933518870d0 &
    +t2*44.1848041560412224d0)) &
    +t*(3479.95663350926514d0 &
    +t2*(482.900172581418890d0 &
    +t2)))
elseif(mc.gt.0.206924d0) then
    t=6.95255575949719117d0*mc-1.43865064797819679d0
    t2=t*t
    ceik=((4870.25402224986382d0 &
    +t2*(3738.29369283392307d0 &
    +t2*51.3609902253065926d0)) &
    +t*(7307.18826377416591d0 &
    +t2*(754.928587580583704d0 &
    +t2*0.571948962277566451d0)))/ &
    ((2172.51745704102287d0 &
    +t2*(2056.13612019430497d0 &
    +t2*44.9026847057686146d0)) &
    +t*(3565.04737778032566d0 &
    +t2*(493.962405117599400d0 &
    +t2)))
elseif(mc.gt.0.121734d0) then
    t=11.7384669562155183d0*mc-1.42897053644793990d0
    t2=t*t
    ceik=((5514.8512729127464d0 &
    +t2*(4313.60788246750934d0 &
    +t2*60.598720224393536d0)) &
    +t*(8350.4595896779631d0 &
    +t2*(880.27903031894216d0 &
    +t2*0.68504458747933773d0)))/ &
    ((2218.41682813309737d0 &
    +t2*(2107.97379949034285d0 &
    +t2*45.6911096775045314d0)) &
    +t*(3650.41829123846319d0 &
    +t2*(505.74295207655096d0 &
    +t2)))
elseif(mc.gt.0.071412d0) then
    t=19.8720241643813839d0*mc-1.41910098962680339d0
    t2=t*t
    ceik=((6188.8743957372448d0 &
    +t2*(4935.41351498551527d0 &
    +t2*70.981049144472361d0)) &
    +t*(9459.3331440432847d0 &
    +t2*(1018.21910476032105d0 &
    +t2*0.81599895108245948d0)))/ &
    ((2260.73112539748448d0 &
    +t2*(2159.68721749761492d0 &
    +t2*46.5298955058476510d0)) &
    +t*(3732.66955095581621d0 &
    +t2*(517.86964191812384d0 &
    +t2)))
elseif(mc.gt.0.041770d0) then
    t=33.7359152553808785d0*mc-1.40914918021725929d0
    t2=t*t
    ceik=((6879.5170681289562d0 &
    +t2*(5594.8381504799829d0 &
    +t2*82.452856129147838d0)) &
    +t*(10615.0836403687221d0 &
    +t2*(1167.26108955935542d0 &
    +t2*0.96592719058503951d0)))/ &
    ((2296.88303450660439d0 &
    +t2*(2208.74949754945558d0 &
    +t2*47.3844470709989137d0)) &
    +t*(3807.37745652028212d0 &
    +t2*(529.79651353072921d0 &
    +t2)))
elseif(mc.gt.0.024360d0) then
    t=57.4382538770821367d0*mc-1.39919586444572085d0
    t2=t*t
    ceik=((7570.6827538712100d0 &
    +t2*(6279.2661370014890d0 &
    +t2*94.886883830605940d0)) &
    +t*(11792.9392624454532d0 &
    +t2*(1325.01058966228180d0 &
    +t2*1.13537029594409690d0)))/ &
    ((2324.04824540459984d0 &
    +t2*(2252.22250562615338d0 &
    +t2*48.2089280211559345d0)) &
    +t*(3869.56755306385732d0 &
    +t2*(540.85752251676412d0 &
    +t2)))
elseif(mc.gt.0.014165d0) then
    t=98.0872976949485042d0*mc-1.38940657184894556d0
    t2=t*t
    ceik=((8247.2601660137746d0 &
    +t2*(6974.7495213178613d0 &
    +t2*108.098282908839979d0)) &
    +t*(12967.7060124572914d0 &
    +t2*(1488.54008220335966d0 &
    +t2*1.32411616748380686d0)))/ &
    ((2340.47337508405427d0 &
    +t2*(2287.70677154700516d0 &
    +t2*48.9575432570382154d0)) &
    +t*(3915.63324533769906d0 &
    +t2*(550.45072377717361d0 &
    +t2)))
elseif(mc.gt.0.008213d0) then
    t=168.010752688172043d0*mc-1.37987231182795699d0
    t2=t*t
    ceik=((8894.2961573611293d0 &
    +t2*(7666.5611739483371d0 &
    +t2*121.863474964652041d0)) &
    +t*(14113.7038749808951d0 &
    +t2*(1654.60731579994159d0 &
    +t2*1.53112170837206117d0)))/ &
    ((2344.88618943372377d0 &
    +t2*(2313.28396270968662d0 &
    +t2*49.5906602613891184d0)) &
    +t*(3942.81065054556536d0 &
    +t2*(558.07615380622169d0 &
    +t2)))
elseif(mc.gt.0.d0) then
    t=1.d0-121.758188238159016d0*mc
    ceik=-log(mc*0.0625d0) &
    *(34813.4518336350547d0 &
    +t*(235.767716637974271d0 &
    +t*0.199792723884069485d0))/ &
    (69483.5736412906324d0 &
    +t*(614.265044703187382d0 &
    +t)) &
    -mc*(9382.53386835986099d0 &
    +t*(51.6478985993381223d0 &
    +t*0.00410754154682816898d0))/ &
    (37327.7262507318317d0 &
    +t*(408.017247271148538d0 &
    +t))
else
    write(*,"(a20,1pe15.7)") "(ceik) out of range: mc=",mc
endif
!
!write(*,"(1x,a20,0p5f20.15)") "(ceik) mc,K=",mc,ceik
!
return
end
!===============================================================================
real*8 function ceie(mc)
!
! Double precision minimax rational approximation of E(m):
! the complete elliptic integral of the second kind
!
!     Reference: T. Fukushima, (2015), J. Comp. Appl. Math., 282, 71-76
!
!      "Precise and Fast Computation of General Complete Elliptic
!       Integral of Second Kind"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
real*8 mc
real*8 t,t2
if(mc.gt.1.d0) then
    write(*,"(a20,1pe15.7)") "(ceie) negative m: mc=",mc
elseif(mc.gt.0.566638d0) then
	t=2.30753965506897236d0*mc-1.30753965506897236d0
	t2=t*t
    ceie=((19702.2363352671642d0 &
    +t2*(18177.1879313824040d0 &
    +t2*409.975559128654710d0)) &
    +t*(31904.1559574281609d0 &
    +t2*(4362.94760768571862d0 &
    +t2*10.3244775335024885d0)))/ &
    ((14241.2135819448616d0 &
    +t2*(10266.4884503526076d0 &
    +t2*117.162100771599098d0)) &
    +t*(20909.9899599927367d0 &
    +t2*(1934.86289070792954d0 &
    +t2)))
elseif(mc.gt.0.315153d0) then
    t=3.97638030101198879d0*mc-1.25316818100483130d0
    t2=t*t
    ceie=((16317.0721393008221d0 &
    +t2*(15129.4009798463159d0 &
    +t2*326.113727011739428d0)) &
    +t*(26627.8852140835023d0 &
    +t2*(3574.15857605556033d0 &
    +t2*7.93163724081373477d0)))/ &
    ((13047.1505096551210d0 &
    +t2*(9964.25173735060361d0 &
    +t2*117.670514069579649d0)) &
    +t*(19753.5762165922376d0 &
    +t2*(1918.72232033637537d0 &
    +t2)))
elseif(mc.gt.0.171355d0) then
    t=6.95419964116329852d0*mc-1.19163687951153702d0
    t2=t*t
    ceie=((13577.3850240991520d0 &
    +t2*(12871.9137872656293d0 &
    +t2*263.964361648520708d0)) &
    +t*(22545.4744699553993d0 &
    +t2*(3000.74575264868572d0 &
    +t2*6.08522443139677663d0)))/ &
    ((11717.3306408059832d0 &
    +t2*(9619.40382323874064d0 &
    +t2*118.690522739531267d0)) &
    +t*(18431.1264424290258d0 &
    +t2*(1904.06010727307491d0 &
    +t2)))
elseif(mc.gt.0.090670d0) then
    t=12.3938774245522712d0*mc-1.12375286608415443d0
    t2=t*t
    ceie=((11307.9485341543712d0 &
    +t2*(11208.6068472959372d0 &
    +t2*219.253495956962613d0)) &
    +t*(19328.6173704569489d0 &
    +t2*(2596.54874477084334d0 &
    +t2*4.66931143174036616d0)))/ &
    ((10307.6837501971393d0 &
    +t2*(9241.7604666150102d0 &
    +t2*120.498555754227847d0)) &
    +t*(16982.2450249024383d0 &
    +t2*(1893.41905403040679d0 &
    +t2)))
elseif(mc.gt.0.046453d0) then
    t=22.6157360291290680d0*mc-1.05056878576113260d0
    t2=t*t
    ceie=((9383.1490856819874d0 &
    +t2*(9977.2498973537718d0 &
    +t2*188.618148076418837d0)) &
    +t*(16718.9730458676860d0 &
    +t2*(2323.49987246555537d0 &
    +t2*3.59313532204509922d0)))/ &
    ((8877.1964704758383d0 &
    +t2*(8840.2771293410661d0 &
    +t2*123.422125687316355d0)) &
    +t*(15450.0537230364062d0 &
    +t2*(1889.13672102820913d0 &
    +t2)))
elseif(mc.gt.0.022912d0) then
    t=42.4790790535661187d0*mc-0.973280659275306911d0
    t2=t*t
    ceie=((7719.1171817802054d0 &
    +t2*(9045.3996063894006d0 &
    +t2*169.386557799782496d0)) &
    +t*(14521.7363804934985d0 &
    +t2*(2149.92068078627829d0 &
    +t2*2.78515570453129137d0)))/ &
    ((7479.7539074698012d0 &
    +t2*(8420.3848818926324d0 &
    +t2*127.802109608726363d0)) &
    +t*(13874.4978011497847d0 &
    +t2*(1892.69753150329759d0 &
    +t2)))
elseif(mc.gt.0.010809d0) then
    t=82.6241427745187144d0*mc-0.893084359249772784d0
    t2=t*t
    ceie=((6261.6095608987273d0 &
    +t2*(8304.3265605809870d0 &
    +t2*159.371262600702237d0)) &
    +t*(12593.0874916293982d0 &
    +t2*(2048.68391263416822d0 &
    +t2*2.18867046462858104d0)))/ &
    ((6156.4532048239501d0 &
    +t2*(7979.7435857665227d0 &
    +t2*133.911640385965187d0)) &
    +t*(12283.8373999680518d0 &
    +t2*(1903.60556312663537d0 &
    +t2)))
elseif(mc.gt.0.004841d0) then
    t=167.560321715817694d0*mc-0.811159517426273458d0
    t2=t*t
    ceie=((4978.06146583586728d0 &
    +t2*(7664.6703673290453d0 &
    +t2*156.689647694892782d0)) &
    +t*(10831.7178150656694d0 &
    +t2*(1995.66437151562090d0 &
    +t2*1.75859085945198570d0)))/ &
    ((4935.56743322938333d0 &
    +t2*(7506.8028283118051d0 &
    +t2*141.854303920116856d0)) &
    +t*(10694.5510113880077d0 &
    +t2*(1918.38517009740321d0 &
    +t2)))
elseif(mc.gt.0.d0) then
    t=1.d0-206.568890725056806d0*mc
    ceie=-mc*log(mc*0.0625d0) &
    *(41566.6612602868736d0 &
    +t*(154.034981522913482d0 &
    +t*0.0618072471798575991d0))/ &
    (165964.442527585615d0 &
    +t*(917.589668642251803d0 &
    +t)) &
    +(132232.803956682877d0 &
    +t*(353.375480007017643d0 &
    -t*1.40105837312528026d0))/ &
    (132393.665743088043d0 &
    +t*(192.112635228732532d0 &
    -t))
else
    write(*,"(a20,1pe15.7)") "(ceie) out of range: mc=",mc
endif
!
!write(*,"(1x,a20,0p5f20.15)") "(ceie) mc,E=",mc,ceie
!
return
end
!===============================================================================
real*8 function ceib(mc)
!
! Double precision minimax rational approximation of B(m):
! the first associate complete elliptic integral of the second kind
!
!     Reference: T. Fukushima, (2015), J. Comp. Appl. Math., 282, 71-76
!
!      "Precise and Fast Computation of General Complete Elliptic
!       Integral of Second Kind"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
real*8 mc
real*8 t,t2
if(mc.gt.1.d0) then
    write(*,"(a20,1pe15.7)") "(ceib) negative m: mc=",mc
elseif(mc.gt.0.555073d0) then
	t=2.24755971204264969d0*mc-1.24755971204264969d0
    t2=t*t
	ceib=((2030.25011505956379d0 &
    +t2*(1727.60635612511943d0 &
    +t2*25.0715510300422010d0)) &
    +t*(3223.16236100954529d0 &
    +t2*(361.164121995173076d0 &
    +t2*0.280355207707726826d0)))/ &
    ((2420.64907902774675d0 &
    +t2*(2327.48464880306840d0 &
    +t2*47.9870997057202318d0)) &
    +t*(4034.28168313496638d0 &
    +t2*(549.234220839203960d0 &
    +t2)))
elseif(mc.gt.0.302367d0) then
    t=3.95716761770595079d0*mc-1.19651690106289522d0
    t2=t*t
    ceib=((2209.26925068374373d0 &
    +t2*(1981.37862223307242d0 &
    +t2*29.7612810087709299d0)) &
    +t*(3606.58475322372526d0 &
    +t2*(422.693774742063054d0 &
    +t2*0.334623999861181980d0)))/ &
    ((2499.57898767250755d0 &
    +t2*(2467.63998386656941d0 &
    +t2*50.0198090806651216d0)) &
    +t*(4236.30953048456334d0 &
    +t2*(581.879599221457589d0 &
    +t2)))
elseif(mc.gt.0.161052d0) then
    t=7.07638962601280827d0*mc-1.13966670204861480d0
    t2=t*t
    ceib=((2359.14823394150129d0 &
    +t2*(2254.30785457761760d0 &
    +t2*35.2259786264917876d0)) &
    +t*(3983.28520266051676d0 &
    +t2*(492.601686517364701d0 &
    +t2*0.396605124984359783d0)))/ &
    ((2563.95563932625156d0 &
    +t2*(2633.23323959119935d0 &
    +t2*52.6711647124832948d0)) &
    +t*(4450.19076667898892d0 &
    +t2*(622.983787815718489d0 &
    +t2)))
elseif(mc.gt.0.083522d0) then
    t=12.8982329420869341d0*mc-1.07728621178898491d0
    t2=t*t
    ceib=((2464.65334987833736d0 &
    +t2*(2541.68516994216007d0 &
    +t2*41.5832527504007778d0)) &
    +t*(4333.38639187691528d0 &
    +t2*(571.53606797524881d0 &
    +t2*0.465975784547025267d0)))/ &
    ((2600.66956117247726d0 &
    +t2*(2823.69445052534842d0 &
    +t2*56.136001230010910d0)) &
    +t*(4661.64381841490914d0 &
    +t2*(674.25435972414302d0 &
    +t2)))
elseif(mc.gt.0.041966d0) then
    t=24.0639137549331023d0*mc-1.00986620463952257d0
    t2=t*t
    ceib=((2509.86724450741259d0 &
    +t2*(2835.27071287535469d0 &
    +t2*48.9701196718008345d0)) &
    +t*(4631.12336462339975d0 &
    +t2*(659.86172161727281d0 &
    +t2*0.54158304771955794d0)))/ &
    ((2594.15983397593723d0 &
    +t2*(3034.20118545214106d0 &
    +t2*60.652838995496991d0)) &
    +t*(4848.17491604384532d0 &
    +t2*(737.15143838356850d0 &
    +t2)))
elseif(mc.gt.0.020313d0) then
    t=46.1829769546944996d0*mc-0.938114810880709371d0
    t2=t*t
    ceib=((2480.58307884128017d0 &
    +t2*(3122.00900554841322d0 &
    +t2*57.541132641218839d0)) &
    +t*(4845.57861173250699d0 &
    +t2*(757.31633816400643d0 &
    +t2*0.62119950515996627d0)))/ &
    ((2528.85218300581396d0 &
    +t2*(3253.86151324157460d0 &
    +t2*66.496093157522450d0)) &
    +t*(4979.31783250484768d0 &
    +t2*(812.40556572486862d0 &
    +t2)))
elseif(mc.gt.0.009408d0) then
    t=91.7010545621274645d0*mc-0.862723521320495186d0
    t2=t*t
    ceib=((2365.25385348859592d0 &
    +t2*(3381.09304915246175d0 &
    +t2*67.442026950538221d0)) &
    +t*(4939.53925884558687d0 &
    +t2*(862.16657576129841d0 &
    +t2*0.70143698925710129d0)))/ &
    ((2390.48737882063755d0 &
    +t2*(3462.34808443022907d0 &
    +t2*73.934680452209164d0)) &
    +t*(5015.4675579215077d0 &
    +t2*(898.99542983710459d0 &
    +t2)))
elseif(mc.gt.0.004136d0) then
    t=189.681335356600910d0*mc-0.784522003034901366d0
    t2=t*t
    ceib=((2160.82916040868119d0 &
    +t2*(3584.53058926175721d0 &
    +t2*78.769178005879162d0)) &
    +t*(4877.14832623847052d0 &
    +t2*(970.53716686804832d0 &
    +t2*0.77797110431753920d0)))/ &
    ((2172.70451405048305d0 &
    +t2*(3630.52345460629336d0 &
    +t2*83.173163222639080d0)) &
    +t*(4916.35263668839769d0 &
    +t2*(993.36676027886685d0 &
    +t2)))
elseif(mc.gt.0.d0) then
    t=1.d0-106.292517006802721d0*mc
    ceib=mc*log(mc*0.0625d0) &
    *(6607.46457640413908d0 &
    +t*(19.0287633783211078d0 &
    -t*0.00625368946932704460d0))/ &
    (26150.3443630974309d0 &
    +t*(354.603981274536040d0 &
    +t)) &
    +(26251.5678902584870d0 &
    +t*(168.788023807915689d0 &
    +t*0.352150236262724288d0))/ &
    (26065.7912239203873d0 &
    +t*(353.916840382280456d0 &
    +t))
else
    write(*,"(a20,1pe15.7)") "(ceib) out of range: mc=",mc
endif
!
!write(*,"(1x,a20,0p5f20.15)") "(ceib) mc,B=",mc,ceib
!
return
end
!===============================================================================
real*8 function ceid(mc)
!
! Double precision minimax rational approximation of D(m):
! the second associate complete elliptic integral of the second kind
!
!     Reference: T. Fukushima, (2015), J. Comp. Appl. Math., 282, 71-76
!
!      "Precise and Fast Computation of General Complete Elliptic
!       Integral of Second Kind"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
real*8 mc
real*8 t,t2
if(mc.gt.1.d0) then
    write(*,"(a20,1pe15.7)") "(ceid) negative m: mc=",mc
elseif(mc.gt.0.599909d0) then
	t=2.49943137936119533d0*mc-1.49943137936119533d0
	t2=t*t
    ceid=((1593.39813781813498d0 &
    +t2*(1058.56241259843217d0 &
    +t2*11.7584241242587571d0)) &
    +t*(2233.25576544961714d0 &
    +t2*(195.247394601357872d0 &
    +t2*0.101486443490307517d0)))/ &
    ((1685.47865546030468d0 &
    +t2*(1604.88100543517015d0 &
    +t2*38.6743012128666717d0)) &
    +t*(2756.20968383181114d0 &
    +t2*(397.504162950935944d0 &
    +t2)))
elseif(mc.gt.0.359180d0) then
    t=4.15404874360795750d0*mc-1.49205122772910617d0
    t2=t*t
    ceid=((1967.01442513777287d0 &
    +t2*(1329.30058268219177d0 &
    +t2*15.0447805948342760d0)) &
    +t*(2779.87604145516343d0 &
    +t2*(247.475085945854673d0 &
    +t2*0.130547566005491628d0)))/ &
    ((1749.70634057327467d0 &
    +t2*(1654.40804288486242d0 &
    +t2*39.1895256017535337d0)) &
    +t*(2853.92630369567765d0 &
    +t2*(406.925098588378587d0 &
    +t2)))
elseif(mc.gt.0.214574d0) then
    t=6.91534237860116454d0*mc-1.48385267554596628d0
    t2=t*t
    ceid=((2409.64196912091452d0 &
    +t2*(1659.30176823041376d0 &
    +t2*19.1942111405094383d0)) &
    +t*(3436.40744503228691d0 &
    +t2*(312.186468430688790d0 &
    +t2*0.167847673021897479d0)))/ &
    ((1824.89205701262525d0 &
    +t2*(1715.38574780156913d0 &
    +t2*39.8798253173462218d0)) &
    +t*(2971.02216287936566d0 &
    +t2*(418.929791715319490d0 &
    +t2)))
elseif(mc.gt.0.127875d0) then
    t=11.5341584101316047d0*mc-1.47493050669557896d0
    t2=t*t
    ceid=((2926.81143179637839d0 &
    +t2*(2056.45624281065334d0 &
    +t2*24.3811986813439843d0)) &
    +t*(4214.52119721241319d0 &
    +t2*(391.420514384925370d0 &
    +t2*0.215574280659075512d0)))/ &
    ((1910.33091918583314d0 &
    +t2*(1787.99942542734799d0 &
    +t2*40.7663012893484449d0)) &
    +t*(3107.04531802441481d0 &
    +t2*(433.673494280825971d0 &
    +t2)))
elseif(mc.gt.0.076007d0) then
    t=19.2797100331611013d0*mc-1.46539292049047582d0
    t2=t*t
    ceid=((3520.63614251102960d0 &
    +t2*(2526.67111759550923d0 &
    +t2*30.7739877519417978d0)) &
    +t*(5121.2842239226937d0 &
    +t2*(486.926821696342529d0 &
    +t2*0.276315678908126399d0)))/ &
    ((2003.81997889501324d0 &
    +t2*(1871.05914195570669d0 &
    +t2*41.8489850490387023d0)) &
    +t*(3259.09205279874214d0 &
    +t2*(451.007555352632053d0 &
    +t2)))
elseif(mc.gt.0.045052d0) then
    t=32.3049588111775157d0*mc-1.45540300436116944d0
    t2=t*t
    ceid=((4188.00087087025347d0 &
    +t2*(3072.05695847158556d0 &
    +t2*38.5070211470790031d0)) &
    +t*(6156.0080960857764d0 &
    +t2*(599.76666155374012d0 &
    +t2*0.352955925261363680d0)))/ &
    ((2101.60113938424690d0 &
    +t2*(1961.76794074710108d0 &
    +t2*43.0997999502743622d0)) &
    +t*(3421.55151253792527d0 &
    +t2*(470.407158843118117d0 &
    +t2)))
elseif(mc.gt.0.026626d0) then
    t=54.2711386084880061d0*mc-1.44502333658960165d0
    t2=t*t
    ceid=((4916.74442376570733d0 &
    +t2*(3688.12811638360551d0 &
    +t2*47.6447145147811350d0)) &
    +t*(7304.6632479558695d0 &
    +t2*(729.75841970840314d0 &
    +t2*0.448422756936257635d0)))/ &
    ((2197.49982676612397d0 &
    +t2*(2055.19657857622715d0 &
    +t2*44.4576261146308645d0)) &
    +t*(3584.94502590860852d0 &
    +t2*(490.880160668822953d0 &
    +t2)))
elseif(mc.gt.0.015689d0) then
    t=91.4327512114839536d0*mc-1.43448843375697175d0
    t2=t*t
    ceid=((5688.7542903989517d0 &
    +t2*(4364.21513060078954d0 &
    +t2*58.159468141567195d0)) &
    +t*(8542.6096475195826d0 &
    +t2*(875.35992968472914d0 &
    +t2*0.56528145509695951d0)))/ &
    ((2285.44062680812883d0 &
    +t2*(2145.80779422696555d0 &
    +t2*45.8427480379028781d0)) &
    +t*(3739.30422133833258d0 &
    +t2*(511.23253971875808d0 &
    +t2)))
elseif(mc.gt.0.009216d0) then
    t=154.487872701992894d0*mc-1.42376023482156651d0
    t2=t*t
    ceid=((6475.3392225234969d0 &
    +t2*(5081.2997108708577d0 &
    +t2*69.910123337464043d0)) &
    +t*(9829.1138694605662d0 &
    +t2*(1033.32687775311981d0 &
    +t2*0.70526087421186325d0)))/ &
    ((2357.74885505777295d0 &
    +t2*(2226.89527217032394d0 &
    +t2*47.1609071069631012d0)) &
    +t*(3872.32565152553360d0 &
    +t2*(530.03943432061149d0 &
    +t2)))
elseif(mc.gt.0.d0) then
    t=1.d0-108.506944444444444d0*mc
    ceid=-log(mc*0.0625d0) &
    *(6.2904323649908115d6 &
    +t*(58565.284164780476d0 &
    +t*(131.176674599188545d0 &
    +t*0.0426826410911220304d0)))/ &
    (1.24937550257219890d7 &
    +t*(203580.534005225410d0 &
    +t*(921.17729845011868d0 &
    +t))) &
    -(27356.1090344387530d0 &
    +t*(107.767403612304371d0 &
    -t*0.0827769227048233593d0))/ &
    (27104.0854889805978d0 &
    +t*(358.708172147752755d0 &
    +t))
else
    write(*,"(a20,1pe15.7)") "(ceid) out of range: mc=",mc
endif
!
!write(*,"(1x,a20,0p5f20.15)") "(ceid) mc,D=",mc,ceid
!
return
end
!===============================================================================
real*8 function ceis(mc)
!
! Double precision minimax rational approximation of S(m):
! the special complete elliptic integral
!
!     Reference: T. Fukushima, (2016), Astron. J., re-revised
!
!      "Zonal Toroidal Harmonic Expansions of External Gravitational
!       Fields for Ring-Like Objects"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
real*8 mc
real*8 t,t2
if(mc.gt.1.d0) then
    write(*,"(a20,1pe15.7)") "(ceis) negative m: mc=",mc
elseif(mc.gt.0.601505d0) then
    t=2.5094417746772230517d0*mc-1.5094417746772230517d0
	t2=t*t
    ceis=((1026.17763535952307d0 &
    +t2*(485.677469611677641d0 &
    -t2*1.03090351230372971d0)) &
    +t*(1274.82870240508851d0 &
    +t2*(52.3124844458365222d0 &
    +t2*0.00739324234515133022d0)))/ &
    ((3634.91161208781745d0 &
    +t2*(3687.92980597191086d0 &
    +t2*72.4741783934871057d0)) &
    +t*(6150.43453804219490d0 &
    +t2*(908.926554921646573d0 &
    -t2)))
elseif(mc.gt.0.399058d0) then
    t=4.9395644292086323828d0*mc-1.9711727019911384214d0
    t2=t*t
    ceis=((86.7180733051785905d0 &
    +t2*(-36.6337761424801254d0 &
    -t2*1.38891143351659923d0)) &
    +t*(15.7273307842061126d0 &
    +t2*(-14.8291225844217536d0 &
    +t2*0.00361918439375104252d0)))/ &
    ((235.196742248774012d0 &
    +t2*(-88.9257271644547379d0 &
    -t2*15.3339604566201331d0)) &
    +t*(116.962483360689959d0 &
    +t2*(-71.2170053927684721d0 &
    -t2)))
elseif(mc.gt.0.219501d0) then
    t=5.569262128460600255d0*mc-1.2224586064592302166d0
    t2=t*t
    ceis=((74.7877514681889701d0 &
    +t2*(82.9949360943443884d0 &
    +t2*1.56773125483764527d0)) &
    +t*(134.483076192401992d0 &
    +t2*(20.2462516310088192d0 &
    -t2*0.00367147603407009175d0)))/ &
    ((143.112254337953479d0 &
    +t2*(267.680587804708068d0 &
    +t2*17.1824818506270368d0)) &
    +t*(321.171975806674627d0 &
    +t2*(101.689772698952843d0 &
    +t2)))
elseif(mc.gt.0.128079d0) then
    t=11.059255490920351242d0*mc-1.4275176395125080180d0
    t2=t*t
    ceis=((368.838713207397348d0 &
    +t2*(273.591392228690348d0 &
    +t2*3.11283723054624760d0)) &
    +t*(546.423411840019851d0 &
    +t2*(52.8283480246527158d0 &
    -t2*0.00267666613712733130d0)))/ &
    ((537.483931788056883d0 &
    +t2*(649.569549167506858d0 &
    +t2*24.2918898835263356d0)) &
    +t*(977.761239499760638d0 &
    +t2*(191.900997942387173d0 &
    +t2)))
elseif(mc.gt.0.071412d0) then
    t=19.090527280363483639d0*mc-1.4641861708220381047d0
    t2=t*t
    ceis=((737.51190337498540d0 &
    +t2*(509.23908689777924d0 &
    +t2*5.2334058551069839d0)) &
    +t*(1056.48245489063616d0 &
    +t2*(93.951657090929254d0 &
    +t2*0.00127302402727680993d0)))/ &
    ((848.68441046566805d0 &
    +t2*(911.51082544227172d0 &
    +t2*28.2934744382382281d0)) &
    +t*(1462.80391382404865d0 &
    +t2*(248.592677538715336d0 &
    +t2)))
elseif(mc.gt.0.045739d0) then
    t=32.301828283480845016d0*mc-1.4774533238581303702d0
    t2=t*t
    ceis=((1228.30248538045050d0 &
    +t2*(830.59901244101972d0 &
    +t2*8.3445760764027038d0)) &
    +t*(1742.14135275252047d0 &
    +t2*(151.363543340919925d0 &
    +t2*0.0109389121873063296d0)))/ &
    ((1148.72064907022826d0 &
    +t2*(1154.29906916591534d0 &
    +t2*31.7183859972304475d0)) &
    +t*(1922.56327226047254d0 &
    +t2*(299.505582648248342d0 &
    +t2)))
elseif(mc.gt.0.027299d0) then
    t=54.22993492407809111d0*mc-1.4804229934924078091d0
    t2=t*t
    ceis=((1850.32817333933852d0 &
    +t2*(1251.96912463087426d0 &
    +t2*12.7282329449538449d0)) &
    +t*(2624.51676940821043d0 &
    +t2*(228.600812362065591d0 &
    +t2*0.0301997195942178328d0)))/ &
    ((1440.47401213691609d0 &
    +t2*(1388.46074898128627d0 &
    +t2*34.8987556516293270d0)) &
    +t*(2368.67335584121873d0 &
    +t2*(347.987916163266950d0 &
    +t2)))
elseif(mc.gt.0.016280d0) then
    t=90.75233687267447137d0*mc-1.4774480442871403939d0
    t2=t*t
    ceis=((2594.88367393128139d0 &
    +t2*(1777.35557927706750d0 &
    +t2*18.6194756219725130d0)) &
    +t*(3700.86979543655025d0 &
    +t2*(327.656946454522957d0 &
    +t2*0.063865181191751846d0)))/ &
    ((1716.03161920071730d0 &
    +t2*(1612.19856579654867d0 &
    +t2*37.9050752099373943d0)) &
    +t*(2793.09781104785086d0 &
    +t2*(394.278792974606013d0 &
    +t2)))
elseif(mc.gt.0.009692d0) then
    t=151.79113539769277474d0*mc-1.4711596842744383728d0
    t2=t*t
    ceis=((3441.67541433275474d0 &
    +t2*(2401.49730335565496d0 &
    +t2*26.1675237391586565d0)) &
    +t*(4950.25801126609305d0 &
    +t2*(448.806533020458946d0 &
    +t2*0.117111370758512406d0)))/ &
    ((1967.03500815730011d0 &
    +t2*(1820.81267345256858d0 &
    +t2*40.7252758168263498d0)) &
    +t*(3184.72957744623299d0 &
    +t2*(437.723184552594083d0 &
    +t2)))
elseif(mc.gt.0.d0) then
    t=1.d0-103.17787866281469253d0*mc
    ceis=-log(mc*0.0625d0) &
    *(1.94180589378811075d6 &
    +t*(647.22748042899378d0 &
    +t*(-10.1649835519921633d0 &
    +t*0.00183293824579460819d0)))/ &
    (3.79948492107705041d6 &
    +t*(84831.477117968440d0 &
    +t*(568.51804083313882d0 &
    +t))) &
    -(21159.8549313831840d0 &
    +t*(26.4758907608130305d0 &
    +t*0.00284022078125062437d0))/ &
    (10388.4217213415863d0 &
    +t*(203.745109840801902d0 &
    +t))
else
    write(*,"(a20,1pe15.7)") "(ceis) out of range: mc=",mc
endif
!
!write(*,"(1x,a20,0p5f20.15)") "(ceis) mc,K=",mc,ceis
!
return
end
!===============================================================================
         m                     K(m)                     E(m)                     B(m)                     D(m)                     S(m)
   0.95000    2.908337248444552E+00    1.060473727766278E+00    9.632177529937374E-01    1.945119495450815E+00    1.033580781533765E+00
   0.90000    2.578092113348173E+00    1.104774732704073E+00    9.410728015213955E-01    1.637019311826777E+00    7.732739003393133E-01
   0.85000    2.389016486325580E+00    1.143395791883166E+00    9.235803752168575E-01    1.465436111108722E+00    6.374773363433707E-01
   0.80000    2.257205326820854E+00    1.178489924327839E+00    9.088110737045846E-01    1.348394253116269E+00    5.494789742646051E-01
   0.75000    2.156515647499643E+00    1.211056027568460E+00    8.959028209247316E-01    1.260612826574912E+00    4.862800075335734E-01
   0.70000    2.075363135292470E+00    1.241670567945823E+00    8.843737533686884E-01    1.190989381923780E+00    4.380223265072740E-01
   0.65000    2.007598398424376E+00    1.270707479650150E+00    8.739200618486431E-01    1.133678336575733E+00    3.996281149647539E-01
   0.60000    1.949567749806026E+00    1.298428035046913E+00    8.643348918741715E-01    1.085232857931854E+00    3.681632767628054E-01
   0.55000    1.898924910271554E+00    1.325024497958230E+00    8.554696151564202E-01    1.043455295115134E+00    3.417921453794791E-01
   0.50000    1.854074677301372E+00    1.350643881047676E+00    8.472130847939789E-01    1.006861592507393E+00    3.192970154268276E-01
   0.45000    1.813883936816983E+00    1.375401971871116E+00    8.394795702706128E-01    9.744043665463695E-01    2.998328806127927E-01
   0.40000    1.777519371491253E+00    1.399392138897432E+00    8.322012900067006E-01    9.453180814845527E-01    2.827919786946301E-01
   0.35000    1.744350597225614E+00    1.422691133490879E+00    8.253235579835159E-01    9.190270392420975E-01    2.677242321673756E-01
   0.30000    1.713889448178791E+00    1.445363064412665E+00    8.188015022917051E-01    8.950879458870858E-01    2.542881453179364E-01
   0.25000    1.685750354812596E+00    1.467462209339427E+00    8.125977729199205E-01    8.731525818926756E-01    2.422192358910202E-01
   0.20000    1.659623598610528E+00    1.489035058095853E+00    8.066808960371524E-01    8.529427025733752E-01    2.313090326811136E-01
   0.15000    1.635256732264580E+00    1.510121832092819E+00    8.010240644528450E-01    8.342326678117351E-01    2.213906890592681E-01
   0.10000    1.612441348720219E+00    1.530757636897763E+00    7.956042304956574E-01    8.168371182245620E-01    2.123288772890451E-01
   0.05000    1.591003453790792E+00    1.550973351780472E+00    7.904014135843953E-01    8.006020402063971E-01    2.040125324400383E-01
   0.00000    1.570796326794897E+00    1.570796326794896E+00    7.853981633974484E-01    7.853981633974485E-01    1.963495408493621E-01
!===============================================================================
