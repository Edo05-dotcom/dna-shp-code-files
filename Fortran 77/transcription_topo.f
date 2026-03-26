program supercoildiffusion
      implicit none
      integer*4 nstep,nDNA,npol,nmax,nequil,nequilrun
      parameter(nstep=3000000,nequilrun=1000000,nDNA=1000,npol=10)
      parameter(nmax=100000,nequil=100000)
      integer*4 i,j,index,imoved,n,ngene,genelength
      integer*4 iup,idwn,iupstream,distance
      integer*4 ndata,ndata2,quantum,trial
      integer*4 seed
      integer*4 ntranscript(nDNA)
      integer*4 engaged(npol)
      integer*4 genetranscribedby(npol)
      integer*4 pol(nDNA),polpos(nDNA)
      integer*4 gene(nDNA),genepos(nDNA),genedir(nDNA)
      double precision transcribed(nmax)
      double precision D,velocity,dt,dx
      double precision pgene
      double precision kon,ka,a
      double precision ktopo
      double precision a1,ran2
      double precision timetranscription,flux0
      double precision totalLk,totalLk2,norm,H
      double precision histo(nDNA)
      double precision timeengaged(npol)
      double precision Lk(nDNA),Lknew(nDNA)
      double precision avLk(nDNA)
      double precision flux(nDNA)
      double precision JoverD
      external ran2

c     parameters
      D=10.d0!1.d0 ! supercoil/linking number diffusion
      dt=0.01d0!0.01d0
      dx=1.d0

      timetranscription=10.d0!100.d0
       
      flux0=1.d0!12.d0!9.d0!6.d0!3.d0!6.d0!1.d0
      kon=0.001d0!0.001d0
      a=100.d0!0.1d0 ! coupling kon/supercoiling

      pgene=0.01d0   ! average density of genes along DNA
      seed=37!612!59!11!-1

      quantum=1000 ! determines how often data dumped

c     initialisation and positioning of genes

      ngene=10
      genelength=30!(2*nDNA)/(3*ngene)
      velocity=3.d0!dfloat(genelength)/timetranscription

      open(unit=2,file='in.dat',status='unknown')
      read(2,*) kon
      read(2,*) flux0
      read(2,*) D
      read(2,*) genelength
      read(2,*) ktopo
      read(2,*) seed
      close(2)

      timetranscription=dfloat(genelength)/velocity
      
      JoverD = flux0 / D
      
c     Position of genes
      do i=1,nDNA
         gene(i)=0
         Lk(i)=0.d0
         avLk(i)=0.d0
      enddo
      do i=1,ngene
         pol(i)=0
         ntranscript(i)=0
      enddo
!      do i=1,ngene
!         index=i*nDNA/ngene-nDNA/(2*ngene)
!         gene(index)=1
!         genepos(i)=index
!         a1=ran2(seed)
!         if(a1.le.0.5d0) then
!            genedir(i)=-1
!         else
!            genedir(i)=1
!         endif
!      enddo
      open(unit=2,file='input/genepos.input',status='unknown')
      do j=1,ngene
         read(2,*) i,genepos(i),genedir(i)
      enddo
      close(2)



      open(unit=11,file='output/genepos.dat',status='replace')
      do i=1,ngene
         write(11,*) i,genepos(i),genedir(i)
c         write(6,*) i,genepos(i),genedir(i)
         call flush(11)
      enddo
      close(11)

      do i=1,nDNA
         flux(i)=0.d0
      enddo

      do j=1,npol
         engaged(j)=0
         timeengaged(j)=0.d0
         genetranscribedby(j)=0
      enddo

      ndata=0
      ndata2=0

      open(unit=10,file='output/polymerasepositions.dat',
     1     status='replace')
      open(unit=20,file='output/supercoilprofile.dat',status='replace')
      open(unit=30,file='output/supercoilvariance.dat',status='replace')
      open(unit=40,file='output/supercoilgenes.dat',status='replace')
      open(unit=45,file='output/averagesupercoilingtranscription.dat',
     1     status='replace')
      open(unit=70,file='output/shannon_vs_JD_topo.dat',
     1         status='replace',position='append')

c     write(6,*) 'ok now start'

c     time stepping
      do n=1,nstep
         
c     evolution polymerases
         do j=1,npol
            if(engaged(j).eq.0) then
               a1=ran2(seed)
               trial=int(a1*ngene)+1!determines which gene is tried 
               if(pol(trial).eq.0) then
                  a1=ran2(seed)
                  i=genepos(trial)
                  iupstream=i-5*genedir(trial)!looks at Lk/supercoil upstream
                  if(iupstream.lt.1) iupstream=iupstream+nDNA
                  if(iupstream.gt.nDNA) iupstream=iupstream-nDNA
                  ka=kon*(1.d0-a*Lk(iupstream))
                  if(a1.le.ka*dt) then
                     write(6,*) 'gene ', trial,' on with Lk ',
     1                    Lk(iupstream), ' and sense ', genedir(trial)
                     if(n.ge.nequilrun) then
                        ndata=ndata+1
                        transcribed(ndata)=trial
                        ntranscript(trial)=ntranscript(trial)+1
                     endif
                     pol(trial)=1
                     genetranscribedby(j)=trial
                     polpos(j)=genepos(trial)
                     engaged(j)=1
                     timeengaged(j)=0.d0
                     flux(polpos(j))=flux0*dfloat(genedir(trial))
                  endif
               endif
            elseif(engaged(j).eq.1) then
               timeengaged(j)=timeengaged(j)+dt
               if(timeengaged(j).lt.timetranscription) then
                  flux(polpos(j))=0.d0
                  polpos(j)=genepos(j)+int(velocity*timeengaged(j))*
     1                 genedir(j)
                  flux(polpos(j))=flux0*dfloat(genedir(j))*
     1                 (1.d0+velocity*timeengaged(j))
               elseif(timeengaged(j).ge.timetranscription) then
                  pol(genetranscribedby(j))=0
                  flux(polpos(j))=0
                  engaged(j)=0
                  polpos(j)=0
               endif
            endif
         enddo

c     evolution supercoiling on DNA
         do i=1,nDNA
            iup=i+1
            idwn=i-1
            if(i.eq.1) idwn=nDNA
            if(i.eq.nDNA) iup=1

            Lknew(i)=Lk(i)+dt*D*(Lk(iup)+Lk(idwn)-2.d0*Lk(i))/dx**2.d0
            Lknew(i)=Lknew(i)+dt/dx*(flux(i)-flux(iup))
            Lknew(i)=Lknew(i)-dt*ktopo*Lk(i)
         enddo

         do i=1,nDNA
            Lk(i)=Lknew(i)
         enddo

c     dump output
         if(mod(n,quantum).eq.0.and.n.ge.nequilrun) then
            ndata2=ndata2+1
            rewind(10)
            do j=1,npol
               write(10,*) j,polpos(j),engaged(j),timeengaged(j)
               call flush(10)
            enddo
            rewind(20)
            totalLk=0.d0
            totalLk2=0.d0
            do i=1,nDNA
               avLk(i)=avLk(i)+Lk(i)
               totalLk=totalLk+Lk(i)
               totalLk2=totalLk2+Lk(i)**2.d0
               write(20,*) i,Lk(i)
               call flush(20)
            enddo
            write(30,*) n,totalLk,totalLk2,
     1           (totalLk2-totalLk**2.d0)/dfloat(nDNA)
            call flush(30)
            rewind(40)
            rewind(45)
            do j=1,ngene               
               i=genepos(j)
               iupstream=i-5*genedir(trial)!looks at Lk/supercoil upstream
               if(iupstream.lt.1) iupstream=iupstream+nDNA
               if(iupstream.gt.nDNA) iupstream=iupstream-nDNA
               write(40,*) j,Lk(iupstream)
               write(45,*) j,ntranscript(j),
     1              avLk(iupstream)/dfloat(ndata2),
     1              avLk(i)/dfloat(ndata2),genelength
               call flush(40)
               call flush(45)
            enddo
            
            
         endif
            
      enddo


      close(10)
      close(20)
      close(30)
      close(40)
      close(45)

      open(unit=50,file='output/history.dat',status='replace')
      open(unit=60,file='output/geneprofile.dat',status='replace')

      do i=1,ngene
         histo(i)=0.d0
      enddo

      do n=1,ndata
         write(50,*) n,transcribed(n),genedir(transcribed(n))
         histo(transcribed(n))=histo(transcribed(n))+1.d0/dfloat(ndata)
         call flush(50)
      enddo
      close(50)

      norm=0.d0
      do i=1,ngene
         write(60,*) i,histo(i),histo(i)-1.d0/dfloat(ngene)
         norm=norm+(histo(i)-1.d0/dfloat(ngene))**2.d0
         call flush(60)
      enddo
      close(60)

      norm=norm/(1.d0-1.d0/dfloat(ngene))
      norm=dsqrt(norm)
      write(6,*) norm

c     Shannon entropy
      H = 0.d0
      do i=1,ngene
         if(histo(i).gt.0.d0) then
            H = H - histo(i)*dlog(histo(i))
         endif
      enddo

      write(6,*) 'Shannon entropy H = ',H

c >>> NEW <<<
c     write entropy vs J/D
      write(70,*) JoverD, H
      call flush(70)
      close(70)

      stop
      end

      FUNCTION ran2(idum)
      implicit double precision (a-h,o-z)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *     NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
                   k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue                              
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
       ran2=min(AM*iy,RNMX)
      return
      END
C     (C) Copr. 1986-92 Numerical Recipes Software ')0.