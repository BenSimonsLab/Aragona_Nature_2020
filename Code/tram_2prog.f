      program tram_2prog
c--------------------------------------------------------------------------
c
c     program to simulate and fit expansion (postnatal or healing) of
c     squamous skin epidermis.
c
c     The model involves a two-progenitor compartment model
c     organized into a square lattice of sites. In steady-state, each
c     site plays host to a renewing cell. Asymmetric division gives
c     rise to a renewing cell and a progenitor, committed to differentiation
c     through one round of terminal division before stratification.
c
c---- Check list of things to be changed for new application
c
c     1. itmax - number of time points (note that some work is required for
c                a zero time point
c     2. changes of format statements - search **change
c     3. parameters - nx,ny,dlrmin,dlrmax,dcymin,dcymax
c     4. readex - entry of data
c      
c--------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c--------------------------------------------------------------------------
c
c     The following parameters are common and need to be updated elsewhere
c
c     lattice size = nx * ny
c     nsupr = number of suprabasal cell spaces
c     nmax = max size of arrays for lattice
c     itmax  = number of time points
c
c--------------------------------------------------------------------------
      parameter (nx=1000,ny=1000,nsupr=4,nmax=nx*ny*(2+nsupr),itmax=4)
c--------------------------------------------------------------------------
c
c     The following parameters are local to this subroutine
c
c     if scan:
c        range of r values runs from dlrmin to dlrmax
c        range of cycle times runs from dcymin to dcymax
c        number of steps through range = nstep
c
c     if run:
c        r = dlrmin
c        cycle time = dcymin
c--------------------------------------------------------------------------
      parameter (dlrmin=0.1d0,dlrmax=0.25d0,dcymin=0.6d0,dcymax=1.2d0,
     +     nstep=10)
c--------------------------------------------------------------------------
c
c     distex = experimental clone size distributions
c     distth = theoretical clone size distributions
c     dav = array containing averages
c     dsq = array containing mean-square differences of theory and experiment
c     dtime = time points
c     nsites = array containing number of sites occupied
c
c--------------------------------------------------------------------------
      real*8 distex(0:nmax,2,0:itmax),distth(0:nmax,2,0:itmax),
     +   dav(5,0:itmax),dsq(nstep,nstep),dtime(0:itmax)
      integer*8 nsites(0:itmax)
c---- read experimental data
      call readex(distex,dtime)
c---- scan or run
      write(*,*) 'scan(0) or run(1)'
      read(*,*) iopt
      if (iopt.gt.0) then
         dlrop=dlrmin
         dcyop=dcymin
         goto 250
      endif
c---- get step size in r and cycle time
      dlrstp=(dlrmax-dlrmin)/dble(nstep-1)
      dcyst=(dcymax-dcymin)/dble(nstep-1)
c---- initialise optimal parameter values
      dlsq=1.0d6
      dlrop=0.0d0
      dcyop=0.0d0
c---- loop over r and cycle values
         do 200 l1=1,nstep
         dlrval=dlrmin+dble(l1-1)*dlrstp
            do 180 l2=1,nstep
            dcycle=dcymin+dble(l2-1)*dcyst
            write(*,'(a2,f6.3,a5,f6.3,a15,f6.3,a1,f6.3)')
     +         'r=',dlrval,' t_c=',dcycle,' optimal r/t_c=',dlrop,'/',
     1         dcyop
c---- input theory for given r and cycle time
            call theory(distth,dav,dlrval,dcycle,nclmax,nclmxt,dtime,
     +         nsites)
c---- calculate mean-square value
            dsq(l1,l2)=0.0d0
               do 160 l3=1,nmax
                  do 140 i1=1,itmax
                  dsq(l1,l2)=dsq(l1,l2)
     +               +(distth(l3,1,i1)-distex(l3,1,i1))**2.0d0
     1               +(distth(l3,2,i1)-distex(l3,2,i1))**2.0d0
 140              continue
 160           continue
               write(*,'(a22,f10.2,a1,f10.2)') 'mean sq. val. run/opt=',
     +            dsq(l1,l2),'/',dlsq
               if (dsq(l1,l2).lt.dlsq) then
                  dlsq=dsq(l1,l2)
                  dlrop=dlrval
                  dcyop=dcycle
               endif
 180        continue
 200     continue
c---- compute theory values at optimal r and cycle time for output
 250     call theory(distth,dav,dlrop,dcyop,nclmax,nclmxt,dtime,nsites)
c-----------------------------------------------------------------------
c     output data files
c-----------------------------------------------------------------------
      open(unit=20,file='dist_basal.dat',status='replace')
         do 510 l1=1,nclmax
c---- basal output  **change for changing itmax
         write(20,'(i5,12f12.6)') l1,
     +      (distex(l1,1,i1),distth(l1,1,i1),i1=1,itmax)
 510     continue
      close(unit=20)
      open(unit=22,file='dist_total.dat',status='replace')
         do 515 l1=1,nclmxt
c---- total output  **change for changing itmax
         write(22,'(i5,12f12.6)') l1,
     +      (distex(l1,2,i1),distth(l1,2,i1),i1=1,itmax)
 515     continue
      close(unit=22)
      open(unit=21,file='average.dat',status='replace')
         do 520 i1=1,itmax
c---- averages
         write(21,'(8f12.6)') dtime(i1),
     +      dav(1,i1),dav(2,i1),dav(3,i1),dav(4,i1),100.0d0*dav(5,i1),
     1      100.0d0*dav(3,i1)*dav(5,i1),dble(nsites(i1))/dble(nsites(1))
 520     continue
      close(unit=21)
      open(unit=19,file='lsq.dat',status='replace')
c---- output least-square distribution
      write(19,'(a19)') 'a=ListContourPlot[{'
         do 800 l1=1,nstep-1
         write(19,'(a3,29(f7.1,a3),f7.5,a3)') '{',
     +      (dsq(l1,l2),',',l2=1,nstep-1),
     1      dsq(l1,nstep),'},'
 800     continue
      write(19,'(a3,29(f7.1,a3),f7.1,a3)') '{',
     +   (dsq(nstep,l2),',',l2=1,nstep-1),
     1   dsq(nstep,nstep),'}},'
      write(19,'(a120)') 'LabelStyle -> {FontFamily -> "Helvetica", 
     +   FontSize -> 14}, PlotLabel -> "least-square statistic",' 
      write(19,'(a120)') 'FrameLabel -> {"division time, 
     +   1/\[Lambda]","r value"},'
      write(19,'(a120)') 'PlotLegends -> BarLegend[All, 20],'
      write(19,'(a31,f8.3,a3,f8.4,a5,f8.3,a3,f8.4,a3)') 
     +   'Contours -> 20, DataRange -> {{',
     1   dcymin,',',dcymax,'},{',dlrmin,',',dlrmax,'}}]'
      close(unit=19)
      stop
      end
c--------------------------------------------------------------------------c
c                                                                          c
c                                                             theory       c
c                                                                          c
c--------------------------------------------------------------------------c
      subroutine theory(distth,dav,dlrval,dcycle,nclmax,nclmxt,dtime,
     +   nsites)
c--------------------------------------------------------------------------
c
c     subroutine to simulate postnatal expansion
c
c     The model involves a two-progenitor compartment model
c     organized into a square lattice of sites. In steady-state, each
c     site plays host to a renewing cell. Asymmetric division gives
c     rise to a renewing cell and a progenitor, committed to differentiation
c     through one round of terminal division before stratification.
c
c--------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c--------------------------------------------------------------------------
c
c     Each crypt comprises 2 basal cells (one renewing and one TA)
c     and nsupr suprabasal cells
c
c     lattice size nx * ny
c     itmax  = number of time points
c     dratio = fraction of cells that are suprabasal
c     dfac   = dratio/(2(1-dratio)) is the probability of suprabasal choice
c              to ensure this ratio
c
c--------------------------------------------------------------------------
      parameter (nx=1000,ny=1000,nsupr=4,nmax=nx*ny*(2+nsupr),itmax=4,
     +   dratio=0.5d0,dfac=dratio/2.0d0/(1.0d0-dratio))
c--------------------------------------------------------------------------
c
c     model parameters
c
c     dexpan = factor expansion of tissue
c     dprob  = 1/dexpan = probability of site occupancy by renewing cells at
c              time zero
c     dcycle = average cell cycle time (units of days)
c     dlrval = fraction of divisions that lead to loss/replacement between 
c              neighbouring units in steady-state
c     dtoff  = time offset (units of days)
c     dind   = basal cell induction probability
c
c--------------------------------------------------------------------------
      parameter (dexpan=2.0d0,dprob=1.0d0/dexpan,dtoff=-3.0d0,
     +   dcy0=4.5d0,dind=1.0d0)
c--------------------------------------------------------------------------
c
c     Model arrays
c
c     ncrypt  = contains clone entries with obvious notation
c     itime   = time points in units of numerical step size
c     dtime   = time points in units of days
c     nclsz   = clone size array (basal, total,renewing)
c     ndistcl = clone size distribtion (basal, total)
c     nclsv   = array to store number of surviving clones at each time point
c
c     dav     = array for average clone sizes
c     dlrv    = array containing potential variable dlrvals
c     dcy     = array containing potential variable cycle times      
c
c--------------------------------------------------------------------------
      integer*8 ncrypt(nx,ny,2+nsupr),itime(-1:itmax),
     +   nclsz(-1:nmax,2,0:itmax),ndistcl(0:nmax,2,0:itmax),
     1   nclsv(0:itmax),nsites(0:itmax)
      real*8 dav(5,0:itmax),dlrv(0:itmax),dcy(0:itmax),dtime(0:itmax),
     +   distth(0:nmax,2,0:itmax),dtime1(0:itmax)
c---- seed random number generator
      call random_seed
c-----------------------------------------------------------------------
c     create clonal mark - start by marking only the self-renewing compartment
c
c     ncrypt=0 cell but no clone
c           >1 clonal index
c           -1 no cell
c-----------------------------------------------------------------------
      nmark=0
         do 50 l1=1,nx
            do 40 l2=1,ny
               do 20 l3=1,2+nsupr
               ncrypt(l1,l2,l3)=0
c---- no clonal mark in progenitor or suprabasal layers
               if (l3.gt.1) goto 20
               call random_number(random)
c---- enter clonal mark in renewing compartment
               if (random.le.dind) then
                  nmark=nmark+1
                  ncrypt(l1,l2,1)=nmark
               endif
 20            continue
 40         continue
 50      continue
c---- assign physical time points (days)
      dtime1(0)=0.0d0
         do 52 i1=1,itmax
         dtime1(i1)=dtime(i1)-dtoff
 52      continue
      if (dtime1(1).le.0.0d0) then
         write(*,*) 'time offset too big ',dtoff
         stop
      endif
c---- enter loss/replacement rate and cycle times values across timepoints
c         do 55 i1=0,itmax
c         dlrv(i1)=dlrval
c         dcy(i1)=dcycle
c 55      continue
c---- enter loss/replacement dr values across timepoints
      dlrv(1)=0.21d0
      dlrv(2)=dlrval
      dlrv(3)=dlrval
      dlrv(4)=dlrval
c      dlrv(2)=0.276d0*dlrval
c      dlrv(3)=0.545d0*dlrval
c      dlrv(4)=0.273d0*dlrval
c      dlrv(5)=0.029d0*dlrval
c      dlrv(6)=0.029d0*dlrval
      dcy(1)=dcy0*dcycle
      dcy(2)=dcy0*dcycle
      dcy(3)=dcy0*dcycle
      dcy(4)=dcy0*dcycle
c      dcy(2)=dcy0*dcycle/1.9d0
c      dcy(3)=dcy0*dcycle/2.6d0
c      dcy(4)=dcy0*dcycle/2.3d0
c      dcy(5)=dcy0*dcycle/1.8d0
c      dcy(6)=dcy0*dcycle/1.3d0
c---- given the cell cycle time, convert times into physical units
      itime(0)=0
         do 70 i1=1,itmax
         dt=0.0d0
            do 60 i2=1,i1
            dt=dt+(dtime1(i2)-dtime1(i2-1))/dcy(i2)
 60         continue
         itime(i1)=int(dt*dble(nx*ny))
 70      continue
c-----------------------------------------------------------------------
c     big loop over time course
c-----------------------------------------------------------------------
         do 400 i1=1,itmax
         write(*,'(a5,f7.3,a4,f7.3,2i10)') 'time=',dtime1(i1),' dr=',
     +      dlrv(i1),itime(i1),itime(i1-1)
         dlr=dlrv(i1)
         nsites(i1)=0
         if (i1.ne.1) goto 108
c-----------------------------------------------------------------------
c     create inflation effect by randomly removing cells
c-----------------------------------------------------------------------
c---- nclst = contains starting number of clones
         nclst=0
            do 106 l1=1,nx
               do 104 l2=1,ny
                  do 102 l3=1,2+nsupr
                  call random_number(random)
                  if (random.gt.dprob) then
c---- no cell present
                     ncrypt(l1,l2,l3)=-1
                  elseif (ncrypt(l1,l2,l3).gt.0) then
                     nclst=nclst+1
                  endif
 102              continue
 104           continue
 106        continue
         if (nclst.eq.0) then
            write(*,*) 'no clones!'
            stop
         endif
         write(*,*) 'starting clone number=',nclst,'out of ',nx*ny,
     +      ' sites'
c---- loop over time points
 108        do 310 nstep=itime(i1-1),itime(i1)
c---- choose random xy coordinate
 110        call random_number(random)
            lx=int(random*nx)+1
            call random_number(random)
            ly=int(random*ny)+1
c---- is a cell actually present?
            if (ncrypt(lx,ly,1).eq.-1) goto 310
c---- intra or inter unit cell loss/replacement
            call random_number(random)
            if (random.le.dlr) goto 200
c-----------------------------------------------------------------------
c     intra unit cell loss/replacement - asymmetric and terminal division
c-----------------------------------------------------------------------
               do 120 n1=1,nsupr
               ncrypt(lx,ly,2+nsupr-n1+1)=ncrypt(lx,ly,2+nsupr-n1)
 120           continue
               do 130 n1=1,nsupr+1
               ncrypt(lx,ly,2+nsupr-n1+1)=ncrypt(lx,ly,2+nsupr-n1)
 130           continue
            goto 310
c-----------------------------------------------------------------------
c     inter unit cell loss/replacement - first get neighbour
c-----------------------------------------------------------------------
 200        call neigh(nx,ny,lx,ly,lx1,ly1)
               do 210 n1=1,nsupr-1
               ncrypt(lx1,ly1,2+nsupr-n1+1)=ncrypt(lx1,ly1,2+nsupr-n1)
 210           continue
               do 220 n1=1,nsupr-1
               ncrypt(lx1,ly1,2+nsupr-n1+1)=ncrypt(lx1,ly1,2+nsupr-n1)
 220           continue
            ncrypt(lx1,ly1,4)=ncrypt(lx1,ly1,1)
            ncrypt(lx1,ly1,3)=ncrypt(lx1,ly1,1)
c---- stem cell loss/replacement
            ncrypt(lx1,ly1,1)=ncrypt(lx,ly,1)
 310        continue
c-----------------------------------------------------------------------
c     obtain clonal distributions - first zero out arrays
c-----------------------------------------------------------------------
            do 320 l1=-1,nmark
            nclsz(l1,1,i1)=0
            nclsz(l1,2,i1)=0
 320        continue
c---- loop over all sites and update clone size
            do 390 lx=1,nx
               do 380 ly=1,ny
                  do 330 nb=1,2
c---- count the number of sites that are occupied
                  if (ncrypt(lx,ly,nb).ne.-1) nsites(i1)=nsites(i1)+1
c---- nclsz(x,1,y) basal size of clone x at time y
                  nclsz(ncrypt(lx,ly,nb),1,i1)=
     +               nclsz(ncrypt(lx,ly,nb),1,i1)+1
 330              continue
c---- nclsz(x,2,y) suprabasal size of clone x at time y
                  do 340 ns=1,nsupr/2
                  call random_number(random)
                  if (random.lt.dfac) then
                     nclsz(ncrypt(lx,ly,2+2*ns-1),2,i1)=
     +                  nclsz(ncrypt(lx,ly,2+2*ns-1),2,i1)+1
                     nclsz(ncrypt(lx,ly,2+2*ns),2,i1)=
     +                  nclsz(ncrypt(lx,ly,2+2*ns),2,i1)+1
                  endif
 340              continue
 380           continue
 390        continue
 400     continue
c-----------------------------------------------------------------------
c     compute clone size distributions
c-----------------------------------------------------------------------
         do 460 i1=0,itmax
c---- zero distribution
            do 440 l1=0,nmax
            ndistcl(l1,1,i1)=0
            ndistcl(l1,2,i1)=0
 440        continue
c---- maximum basal clone size
         nclmax=0
c---- maximum total clone size
         nclmxt=0
c---- maximum renewing clone size
         nclmxr=0
c---- number of surviving clones at each time point
         nclsv(i1)=0
c---- loop over all possible clones
            do 450 l1=1,nmark
            n1=nclsz(l1,1,i1)
c---- abandon clone if no basal cell as with experiment
            if (n1.eq.0) goto 450
            nclsv(i1)=nclsv(i1)+1
            if (n1.gt.nclmax) nclmax=n1
c---- ndistcl(l1,1,i1) - distribution of basal clone size
            ndistcl(n1,1,i1)=ndistcl(n1,1,i1)+1
            n1=nclsz(l1,1,i1)+nclsz(l1,2,i1)
            if (n1.gt.nclmxt) nclmxt=n1
c---- ndistcl(l1,2,i1) - distribution of total clone size
            ndistcl(n1,2,i1)=ndistcl(n1,2,i1)+1
 450     continue
c-----------------------------------------------------------------------
c     compute average clones size including and excluding zeros
c
c     dav(1) = average basal clone size from all clones
c     dav(2) = average total clone size from all clones
c     dav(3) = average basal clone size from surviving clones
c     dav(4) = average total clone size from surviving clones
c     dav(5) = surviving clone fraction
c     
c-----------------------------------------------------------------------
         dav(1,i1)=0.0d0
         dav(2,i1)=0.0d0
            do 455 l1=0,nmax
            distth(l1,1,i1)=
     +         100.0d0*dble(ndistcl(l1,1,i1))/dble(nclsv(i1))
            distth(l1,2,i1)=
     +         100.0d0*dble(ndistcl(l1,2,i1))/dble(nclsv(i1))
c---- average basal clone size
            dav(1,i1)=dav(1,i1)+
     +         dble(l1)*dble(ndistcl(l1,1,i1))/dble(nclst)
c---- average total clone size
            dav(2,i1)=dav(2,i1)+
     +         dble(l1)*dble(ndistcl(l1,2,i1))/dble(nclst)
 455        continue
         dav(3,i1)=dav(1,i1)*dble(nclst)/dble(nclsv(i1))
         dav(4,i1)=dav(2,i1)*dble(nclst)/dble(nclsv(i1))
         dav(5,i1)=dble(nclsv(i1))/dble(nclst)
 460     continue
      return
      end
c--------------------------------------------------------------------------c
c                                                                          c
c                                                             neigh        c
c                                                                          c
c--------------------------------------------------------------------------c
      subroutine neigh(nx,ny,lx,ly,lx1,ly1)
c-----------------------------------------------------------------------
c     find a neighbouring site on a square lattice of dimension nx*ny
c-----------------------------------------------------------------------
      call random_number(random)
      if (random.lt.0.5d0) then
         ly1=ly
         call random_number(random)
         if (random.lt.0.5d0) then
            lx1=mod(lx,nx)+1
         else
            lx1=mod(lx+nx-2,nx)+1
         endif
      else
         lx1=lx
         call random_number(random)
         if (random.lt.0.5d0) then
            ly1=mod(ly,ny)+1
         else
            ly1=mod(ly+ny-2,ny)+1
         endif
      endif
      return
      end
c--------------------------------------------------------------------------c
c                                                                          c
c                                                             readex       c
c                                                                          c
c--------------------------------------------------------------------------c
      subroutine readex(distex,dtime)
c-----------------------------------------------------------------------
c
c     read in experimental data
c
c     distex(l1,1,i1) = basal size distribution at time i1
c     distex(l1,2,i1) = total size distribution at time i1
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nx=1000,ny=1000,nsupr=4,nmax=nx*ny*(2+nsupr),itmax=4)
      real*8 draw(100),distex(0:nmax,2,0:itmax),dtime(0:itmax)
         do 10 l1=0,nmax
            do 5 i1=0,itmax
            distex(l1,1,i1)=0.0d0
            distex(l1,2,i1)=0.0d0
 5          continue
 10      continue
      distex(1,1,0)=100.0d0
      distex(1,2,0)=100.0d0
c---- read in experimental data
      open(unit=18,file='tramdist.dat')
c---- first timings
      read(18,*,end=70) i1,(dtime(it),d1,it=1,itmax)
c---- assume that data does not have a zero time point
      dtime(0)=0.0d0
      i1=1
 20   read(18,*,end=50) (draw(it),it=1,1+2*itmax)
         do 30 l1=1,itmax
         distex(i1,1,l1)=draw(2+(l1-1)*2)
         distex(i1,2,l1)=draw(3+(l1-1)*2)
 30      continue
      i1=i1+1
      goto 20
c---- output raw distribution
 50   write(*,*) 'raw distribution'
      write(*,'(a4,20f14.2)') 'time',(dtime(it),it=1,itmax)
      write(*,'(a4,20a7)') '    ','  basal','  total','  basal',
     +   '  total','  basal','  total','  basal','  total','  basal',
     1   '  total','  basal','  total'
         do 60 l1=1,i1
         write(*,'(i4,20f7.2)')
     +      l1,(distex(l1,1,l2),distex(l1,2,l2),l2=1,itmax)
 60      continue
      close(unit=18)
      return
 70   write(*,*) 'error - no data'
      stop
      end
