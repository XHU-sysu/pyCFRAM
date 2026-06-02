      SUBROUTINE S_R(u0,as,ss,pts,rad_base)

*    Subroutine for the calculation of S-R

        include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
        common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
        common /clouds/ pre(100),plwc(100),pde(100), piwc(100)
        common /rains/ prwc(100)
        common /graups/ pgwc(100)
        common /umcon/ umco2, umch4, umn2o
        common /radiat/ fds(100), fus(100), dts(100),
     1                  fdir(100), fuir(100), dtir(100),
     1                  fd(100), fu(100), ht(100)
        common /radnew/ fdsljh(100), fusljh(100),
     1                  fdirljh(100), fuirljh(100),
     1                  fuljh(100),fdljh(100),
     1                  fsljh(100),firljh(100)

        dimension as(mbs), ee(mbir)
        real  rad_base(nv1)
        data ee / mbir * 1.0 /

        call rad (as,u0,ss,pts,ee)

c        print*,fds(1)-fus(1)-fuir(1),fus(1),fuir(1),"clear sky"

        do l=1,nv
           rad_base(l)=ht(l)
        enddo
        rad_base(nv1)=fd(nv1)-fu(nv1)

       return
       end


      SUBROUTINE S_R_MC_cloud(u0,as,ss,pts,rad_base,area_c,sw_base,
     &                      lw_base)

*    Subroutine for the calculation of S-R

        include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
        common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
        common /clouds/ pre(100),plwc(100),pde(100), piwc(100)
        common /rains/ prwc(100)
        common /graups/ pgwc(100)
        common /umcon/ umco2, umch4, umn2o
        common /radiat/ fds(100), fus(100), dts(100),
     1                  fdir(100), fuir(100), dtir(100),
     1                  fd(100), fu(100), ht(100)
        common /radnew/ fdsljh(100), fusljh(100),
     1                  fdirljh(100), fuirljh(100),
     1                  fuljh(100),fdljh(100),
     1                  fsljh(100),firljh(100)

      common /clearsky/fuir_clr(100),fdir_clr(100),
     1                 fus_clr(100),fds_clr(100)
      common /totalsky/fuir_tot(100),fdir_tot(100),
     1                 fus_tot(100),fds_tot(100)

      dimension as(mbs), ee(mbir), area_c(100)
      real  rad_base(nv1), fuir_conv(nv1),fdir_conv(nv1)
      real                 fus_conv(nv1),fds_conv(nv1)
      real sw_base(nv1),lw_base(nv1)
      data ee / mbir * 1.0 /

      call rad (as,u0,ss,pts,ee)

      do l = 1, nv
         clr=fuir_clr(l+1)-fuir_clr(l)
         cld=fuir(l+1)-fuir(l)
         fuir_conv(l)=(1.-area_c(l))*clr+area_c(l)*cld

         clr=fdir_clr(l)-fdir_clr(l+1)
         cld=fdir(l)-fdir(l+1)
         fdir_conv(l)=(1.-area_c(l))*clr+area_c(l)*cld

         clr=fus_clr(l+1)-fus_clr(l)
         cld=fus(l+1)-fus(l)
         fus_conv(l)=(1.-area_c(l))*clr+area_c(l)*cld

         clr=fds_clr(l)-fds_clr(l+1)
         cld=fds(l)-fds(l+1)
         fds_conv(l)=(1.-area_c(l))*clr+area_c(l)*cld
      enddo

      fds_tot(1)=fds(1)    !incoming solar radiation which is not influenced by clouds
      fdir_tot(1)=fdir(1)  ! should be zero at all time
      do l = 1, nv
         fds_tot(l+1)=fds_tot(l)-fds_conv(l)
         fdir_tot(l+1)=fdir_tot(l)-fdir_conv(l)
      enddo

      fuir_tot(nv1)=fuir(nv1)   !upward LW wave radiation at surface is not influenced by clouds
      fus_tot(nv1) = as(1)*fds_tot(nv1)    !upward SW radiation at surface = reflected solar rad.

      do l = nv, 1, -1
         fuir_tot(l)=fuir_tot(l+1)-fuir_conv(l)
         fus_tot(l)=fus_tot(l+1)-fus_conv(l)
      enddo

      do l=1,nv
         rad_base(l) = fds_conv(l)+fus_conv(l)+fdir_conv(l)+fuir_conv(l)
         sw_base(l)=fds_conv(l)+fus_conv(l)
         lw_base(l)=fdir_conv(l)+fuir_conv(l)
      enddo
      l = nv1
      rad_base(l)=fds_tot(l)-fus_tot(l)+fdir_tot(l)-fuir_tot(l)
      sw_base(l)=fds_tot(l)-fus_tot(l)
      lw_base(l)=fdir_tot(l)-fuir_tot(l)
      return
      end


      SUBROUTINE S_R_cloudy(u0,as,ss,pts,rad_base,area_c,sw_base,
     &                      lw_base,water_c,ice_c,iseed,no_cloud_out)

*    Subroutine for the calculation of S-R

        include 'para.file'
        real water_c(100),ice_c(100)

        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
        common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
        common /clouds/ pre(100),plwc(100),pde(100), piwc(100)
        common /rains/ prwc(100)
        common /graups/ pgwc(100)
        common /umcon/ umco2, umch4, umn2o
        common /radiat/ fds(100), fus(100), dts(100),
     1                  fdir(100), fuir(100), dtir(100),
     1                  fd(100), fu(100), ht(100)
        common /radnew/ fdsljh(100), fusljh(100),
     1                  fdirljh(100), fuirljh(100),
     1                  fuljh(100),fdljh(100),
     1                  fsljh(100),firljh(100)

      common /clearsky/fuir_clr(100),fdir_clr(100),
     1                 fus_clr(100),fds_clr(100)
      common /totalsky/fuir_tot(100),fdir_tot(100),
     1                 fus_tot(100),fds_tot(100)

      dimension as(mbs), ee(mbir), area_c(100)

      integer no_cloud(100,100),iseed,iarea_c(nv1),mran
      integer yes_c(nv1),no_c(nv1),no_cloud_out(100,100)

      real  rad_base(nv1)
      real sw_base(nv1),lw_base(nv1)
      data ee / mbir * 1.0 /

c     Apple-to-apple OLD CFRAM raw/CFRAM.zip: GW-base.f generates the
c     no_cloud pattern via ran3 and writes it to disk (unit 52); GW-co2.f,
c     GW-wv.f, GW-o3.f, GW-cloud.f, GW-solar.f, GW-albedo.f, GW-ts.f,
c     GW-warm.f, GW-drdt.f all read this saved pattern and copy it into
c     no_cloud (their S_R_cloudy variants have the ran3 generation block
c     commented out). This makes ALL partials and drdt use the IDENTICAL
c     base cloud-overlap pattern, eliminating MC variance from frc_X.
c
c     We mirror that semantics with a /cloud_fix/ common block:
c       use_fixed_cloud = 0 → generate via ran3 (case 0 only)
c       use_fixed_cloud = 1 → copy from no_cloud_fixed (all subsequent cases)
c     The driver (cfram_fu_1col) sets the flag and saves no_cloud_fixed
c     after case 0.
      common /cloud_fix/ use_fixed_cloud, no_cloud_fixed
      integer use_fixed_cloud, no_cloud_fixed(100,100)

      do l = 1, nv
         iarea_c(l)=int(area_c(l)*100)
c         print*,l,iarea_c(l),area_c(l)
      enddo

      if (use_fixed_cloud .ne. 0) then
c        Use saved base pattern (apple-to-apple OLD CFRAM partial cases)
         no_cloud(:,:) = no_cloud_fixed(:,:)
         no_cloud_out(:,:) = no_cloud_fixed(:,:)
      else
c        Generate via ran3 (apple-to-apple OLD GW-base.f)
         yes_c(:)=0
         no_c(:)=0
         no_cloud(:,:)=99
         no_cloud_out(:,:)=-999

         do igrid = 1, 100
            do l = 1, nv
               if (iarea_c(l) .eq. 0) then
                  no_c(l)=no_c(l)+1
                  no_cloud(igrid,l)=1
                  go to 10
               endif
               mran=int(100*ran3(iseed))
               if(mran .lt. iarea_c(l))then
                  if (yes_c(l) .lt. iarea_c(l))then
                     no_cloud(igrid,l)=0
                     yes_c(l) = yes_c(l) + 1
                  else
                     no_c(l)=no_c(l)+1
                     no_cloud(igrid,l)=1  !regardless of mran, no clouds
                  endif
               else
                  if (no_c(l) .lt. (100 - iarea_c(l)))then
                     no_cloud(igrid,l)=1
                     no_c(l)=no_c(l)+1
                  else
                     yes_c(l)=yes_c(l)+1
                     no_cloud(igrid,l)=0
                  endif
               endif
 10            continue
            enddo
         enddo

         no_cloud_out(:,:)=no_cloud(:,:)
      endif

c       print*,"no_cloud",no_cloud(1,:)
c       print*,"no_cloud_out",no_cloud_out(1,:)
c      do l = 1, nv
c         print*,yes_c(l),yes_c(l)+no_c(l),iarea_c(l),l
c         test=0.0
c         do igrid = 1, 100
c            test=test+no_cloud(igrid,l)
c         enddo
c         print*,test,no_c(l)
c      enddo

c      do l1= 1, 5
c         l2 = (l1-1)*20+1
c         l3 = l2 + 19
c         print*,l2,l3
c         do l = 1, nv
c            write(6,*)"l=",l,(no_cloud(i,l),i=l2, l3)
c         enddo
c      enddo

      rad_base(:)=0.0
      sw_base(:)=0.0
      lw_base(:)=0.0
      fus_tot(:)=0.0
      fds_tot(:)=0.0
      fuir_tot(:)=0.0
      fdir_tot(:)=0.0

      do igrid = 1, 100
         test=0.0
         do l = 1, nv
            test=test+no_cloud(igrid,l)
            if(no_cloud(igrid,l) .eq. 0)then
             pre(l) = 15  !10-15 micron meters for water clouds (assumed)
             pde(l) = 40  !40-45 micro meters for ice clouds
c             pmean=0.5*(pp(l)+pp(l+1))
             pmean = pp(l+1)
             dry_density = pmean*100./(287*pt(l))  ! convert back to pascal
             plwc(l) = 1000.*water_c(l)*dry_density  !converting to g/m^3 from kg/kg
             piwc(l) = 1000.*ice_c(l)*dry_density
             pgwc(l)=0.0
c               print*,l,no_cloud(igrid,l),igrid
            else
               pre(l) = 0.0
               plwc(l) = 0.0
               pde(l) = 0.0
               piwc(l) = 0.0
               prwc(l) = 0.0
               pgwc(l) = 0.0
c               print*,l,no_cloud(igrid,l),"no clouds"
            endif
         enddo

         call rad (as,u0,ss,pts,ee)

c         print*,fds(1)-fus(1)-fuir(1),fus(1),fuir(1),test,igrid
         do l=1,nv
            rad_base(l)=rad_base(l)+ht(l)/100.0
            sw_base(l)=sw_base(l)+fsljh(l)/100.0
            lw_base(l)=lw_base(l)+firljh(l)/100.0
            fus_tot(l)=fus_tot(l)+fus(l)/100.0
            fds_tot(l)=fds_tot(l)+fds(l)/100.0
            fdir_tot(l)=fdir_tot(l)+fdir(l)/100.
            fuir_tot(l)=fuir_tot(l)+fuir(l)/100.
         enddo
         rad_base(nv1)=rad_base(nv1)+
     &                     (fd(nv1)-fu(nv1))/100.0
         sw_base(nv1)=sw_base(nv1)+fsljh(nv1)/100.0
         lw_base(nv1)=lw_base(nv1)+firljh(nv1)/100.0
         l = nv1
         fus_tot(l)=fus_tot(l)+fus(l)/100.0
         fds_tot(l)=fds_tot(l)+fds(l)/100.0
         fdir_tot(l)=fdir_tot(l)+fdir(l)/100.
         fuir_tot(l)=fuir_tot(l)+fuir(l)/100.
      enddo

c      print*,fds_tot(1)-fus_tot(1)-fuir_tot(1),fus_tot(1),fuir_tot(1)
      return
      end

      function ran3(idum)  !uniform random generator between (0 and 1)
      integer idum         !set idum = any negative value to initialize the sequence
      integer mbig,mseed,mz
      real ran3, fac
      parameter(mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
      integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(55)
      save iff,inext,inextp,ma
      data iff /0/
      if(idum .lt. 0 .or. iff .eq. 0)then !initialization
         iff=1
c        mj = mseed - iabs(int(real(idum)))  !initialize mas(55) using the seed idum and the large number mseed
         mj = mseed - iabs(idum) !initialize mas(55) using the seed idum and the large number mseed
         mj = mod(mj,mbig)
         ma(55)=mj
         mk=1
         do 11 i = 1, 54
            ii=mod(21*i,55) !now intialize the rest of the table, in a slightly random order
            ma(ii)=mk       !with numbers that are not especially random
            mk=mj-mk
            if(mk .lt. mz)mk = mk + mbig
            mj = ma(ii)
 11      continue
         do 13 k = 1, 4  !we randomize by "warming up the generator"
            do 12 i = 1, 55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
 12         continue
 13         continue
         inext=0 !prepare indices for our first generated number
         inextp=31  !31 is special and cannot be changed.
         idum=1
      endif
      inext=inext+1 !here is where we start, except on initialization.  increment inext, wrapping around 56 to 1
      if(inext .eq. 56)inext = 1
      inextp=inextp+1   !ditto for inextp
      if(inextp .eq. 56)inextp = 1
      mj=ma(inext)-ma(inextp) !now generate a new random number subtractively.
      if(mj .lt. mz)mj = mj+mbig  !be sure that it is in range
      ma(inext)=mj                 !store it
      ran3=mj*fac                   !and output the derived uniform deviate.
      return
      end
