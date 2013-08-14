module mk
  implicit none
  private
  public :: mk_ass
  integer :: n
  double precision, dimension(:), allocatable :: lbworker, lbjob, msvbj
  double precision, dimension(:,:), allocatable :: cm
  integer, dimension(:), allocatable :: mjbw, mwbj, pwbcj, mswbj
  logical, dimension(:), allocatable :: comworkers
  
  contains
  
  subroutine vec_sub_min(vec)
    double precision, dimension(n) :: vec
    integer :: i
    double precision :: tmp
    
    tmp = vec(1)
    do i = 2, n
      if (vec(i) .ge. tmp) cycle
      tmp = vec(i)
    end do
    vec(:) = vec(:) - tmp
  end subroutine
  
  subroutine reduce()
    integer :: i
    
    do i = 1, n
      call vec_sub_min(cm(i,:))
    end do
    
    do i = 1, n
      call vec_sub_min(cm(:,i))
    end do
  end subroutine
  
  !~subroutine init_feas_sol()
  !~  integer :: i, j
  !~  lbjob(:) = 999.d0
  !~  do i = 1, n
  !~    do j = 1, n
  !~      if (cm(i,j) .ge. lbjob(j)) cycle
  !~      lbjob(j) = cm(i,j)
  !~    end do
  !~  end do
  !~end subroutine
  !~
  subroutine greedy_match()
    integer :: w, j
    do w = 1, n
      do j = 1, n
        if (mjbw(w) .eq. -1 .and. mwbj(j) .eq. -1 .and. cm(w,j) - lbworker(w) - lbjob(j) .eq. 0.d0) then
          call match(w,j)
        end if
      end do
    end do
  end subroutine
  
  subroutine match(w,j)
    integer :: w, j
    mjbw(w) = j
    mwbj(j) = w
  end subroutine
  
  function fetch_um_worker() result(w)
    integer :: w
    do w = 1, n
      if (mjbw(w) .eq. -1) return
    end do
    w = 0
  end function
  
  subroutine init_phase(w)
    integer :: w
    integer :: j
    
    comworkers(:) = .false.
    pwbcj(:) = -1
    comworkers(w)= .true.
    do j = 1, n
      msvbj(j) = cm(w,j) - lbworker(w) - lbjob(j)
      mswbj(j) = w
    end do
  end subroutine
  
  subroutine exec_phase()
    integer :: j
    integer :: msworker, msjob, cjob, pworker, worker, tmp
    double precision :: msvalue, slack
    
    do
      msworker = -1
      msjob = -1
      msvalue = 1.d11
      
      do j = 1, n
        if (pwbcj(j) .eq. -1) then
          if (msvbj(j) .lt. msvalue) then
            msvalue = msvbj(j)
            msworker = mswbj(j)
            msjob = j
          end if
        end if
      end do
        
      if (msvalue .gt. 0.d0) then
        call update_label(msvalue)
      end if
      
      pwbcj(msjob) = msworker
      
      if (mwbj(msjob) .eq. -1) then
        cjob = msjob
        pworker = pwbcj(cjob)
        
        do
          tmp = mjbw(pworker)
          call match(pworker,cjob)
          cjob = tmp
          
          if (cjob .eq. -1) exit
          
          pworker = pwbcj(cjob)
        end do
        
        return
      else
        worker = mwbj(msjob)
        comworkers(worker) = .true.
        
        do j = 1, n
          if (pwbcj(j) .eq. -1) then
            slack = cm(worker,j) - lbworker(worker) - lbjob(j)
            
            if (msvbj(j) .gt. slack) then
              msvbj(j) = slack
              mswbj(j) = worker
            end if
          end if
        end do
      end if
      
    end do
  end subroutine
  
  subroutine update_label(slack)
    double precision :: slack
    integer :: w, j
    
    do w = 1, n
      if (comworkers(w)) then
        lbworker(w) = lbworker(w) + slack
      end if
    end do
    
    do j = 1, n
      if (pwbcj(j) .ne. -1) then
        lbjob(j) = lbjob(j) - slack
      else
        msvbj(j) = msvbj(j) - slack
      end if
    end do
  end subroutine

  subroutine mk_ass(nin,mat,sol)
    integer :: nin
    integer, dimension(nin) :: sol
    double precision, dimension(nin,nin) :: mat
    
    integer :: um_w
    
    n = nin
    allocate(lbworker(n),lbjob(n),msvbj(n))
    allocate(cm(n,n))
    allocate(mjbw(n),mwbj(n),pwbcj(n),mswbj(n))
    allocate(comworkers(n))
    cm = mat
    lbworker(:) = 0.d0
    lbjob(:) = 0.d0
    mjbw(:) = -1
    mwbj(:) = -1
    mswbj(:) = 0
    
    call reduce()
    call greedy_match()
    
    do
      um_w = fetch_um_worker()
      
      if (um_w .eq. 0) exit
      
      call init_phase(um_w)
      call exec_phase()
    end do
    
    sol = mjbw
    
    deallocate(lbworker,lbjob,msvbj)
    deallocate(cm)
    deallocate(mjbw,mwbj,pwbcj,mswbj)
    deallocate(comworkers)
  end subroutine
  
end module

