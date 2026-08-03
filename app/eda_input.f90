! This file is part of gfnff.
! SPDX-License-Identifier: LGPL-3.0-or-later
!
! Minimal EDA-GFNFF input reader.
! The Gaussian section parser is a small Fortran port of the state-machine and
! atom-metadata ideas used by banelib's INPUT/GJF API, supplied by its author.
module eda_gfnff_input
  use iso_fortran_env, only : wp => real64
  implicit none
  private

  real(wp),parameter :: angstrom_to_bohr = 1.0_wp/0.529177249_wp
  integer,parameter :: line_len = 4096

  character(len=2),parameter :: pse(118) = [ &
   & 'H ',                                                                                'He', &
   & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
   & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
   & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
   & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
   & 'Cs','Ba','La',                                                                            &
   &                'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',      &
   &                'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
   & 'Fr','Ra','Ac',                                                                            &
   &                'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',      &
   &                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' ]

  type,public :: eda_input_data
    integer :: nat = 0
    integer :: charge = 0
    integer :: spin = 1
    integer :: nfrag = 0
    integer,allocatable :: at(:)
    integer,allocatable :: fragment(:)
    integer,allocatable :: fragment_charge(:)
    real(wp),allocatable :: xyz(:,:)
    character(len=8) :: format = ''
  end type eda_input_data

  public :: read_eda_input,atomic_symbol

contains

  pure function atomic_symbol(z) result(symbol)
    integer,intent(in) :: z
    character(len=2) :: symbol

    if (z >= 1 .and. z <= size(pse)) then
      symbol = pse(z)
    else
      symbol = 'X '
    end if
  end function atomic_symbol

  subroutine read_eda_input(filename,fragment_spec,charge_spec,data,error_message)
    character(len=*),intent(in) :: filename
    character(len=*),intent(in) :: fragment_spec
    character(len=*),intent(in) :: charge_spec
    type(eda_input_data),intent(out) :: data
    character(len=:),allocatable,intent(out) :: error_message
    character(len=:),allocatable :: ext

    error_message = ''
    ext = file_extension(filename)
    select case (ext)
    case ('gjf','com')
      call read_gaussian(filename,data,error_message)
    case ('xyz')
      call read_xyz(filename,fragment_spec,charge_spec,data,error_message)
    case default
      error_message = 'unsupported input extension: .'//ext//' (expected .xyz, .gjf, or .com)'
    end select
  end subroutine read_eda_input

  subroutine read_xyz(filename,fragment_spec,charge_spec,data,error_message)
    character(len=*),intent(in) :: filename,fragment_spec,charge_spec
    type(eda_input_data),intent(out) :: data
    character(len=:),allocatable,intent(out) :: error_message
    integer :: unit,io,i,z,frag_value
    character(len=line_len) :: line,comment
    character(len=256) :: tokens(32)
    integer :: ntok
    logical :: has_inline,all_inline,ok
    integer,allocatable :: supplied_frag(:),charges(:)
    real(wp) :: x,y,zcoord

    error_message = ''
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=io)
    if (io /= 0) then
      error_message = 'cannot open XYZ file: '//trim(filename)
      return
    end if
    read(unit,*,iostat=io) data%nat
    if (io /= 0 .or. data%nat < 1) then
      close(unit)
      error_message = 'invalid XYZ atom count'
      return
    end if
    read(unit,'(a)',iostat=io) comment
    if (io /= 0) then
      close(unit)
      error_message = 'missing XYZ comment line'
      return
    end if
    call parse_xyz_comment(comment,data%charge,data%spin)

    allocate(data%at(data%nat),data%xyz(3,data%nat),data%fragment(data%nat))
    data%fragment = 0
    all_inline = .true.
    do i=1,data%nat
      read(unit,'(a)',iostat=io) line
      if (io /= 0) then
        close(unit)
        error_message = 'XYZ ended before all atoms were read'
        return
      end if
      call whitespace_tokens(line,tokens,ntok)
      if (ntok < 4) then
        close(unit)
        error_message = 'invalid XYZ atom line '//integer_string(i)
        return
      end if
      z = element_to_z(tokens(1))
      if (z < 1) then
        close(unit)
        error_message = 'unknown element on XYZ atom line '//integer_string(i)//': '//trim(tokens(1))
        return
      end if
      read(tokens(2),*,iostat=io) x
      if (io == 0) read(tokens(3),*,iostat=io) y
      if (io == 0) read(tokens(4),*,iostat=io) zcoord
      if (io /= 0) then
        close(unit)
        error_message = 'non-numeric coordinate on XYZ atom line '//integer_string(i)
        return
      end if
      data%at(i)=z
      data%xyz(:,i)=[x,y,zcoord]*angstrom_to_bohr
      has_inline = .false.
      if (ntok >= 5) then
        call parse_fragment_token(tokens(5),frag_value,has_inline)
      end if
      if (has_inline) then
        data%fragment(i)=frag_value
      else
        all_inline=.false.
      end if
    end do
    close(unit)

    if (len_trim(fragment_spec) > 0) then
      call parse_integer_spec(fragment_spec,data%nat,supplied_frag,ok,error_message)
      if (.not.ok) return
      if (size(supplied_frag) /= data%nat) then
        error_message = '--frag must contain exactly one fragment id per atom'
        return
      end if
      data%fragment = supplied_frag
    else if (.not.all_inline) then
      error_message = 'XYZ requires fragment ids in column 5 or via --frag'
      return
    end if

    call normalize_fragments(data%fragment,data%nfrag,error_message)
    if (len(error_message) > 0) return
    if (data%nfrag < 2) then
      error_message = 'EDA-GFNFF requires at least two fragments'
      return
    end if

    if (len_trim(charge_spec) == 0) then
      error_message = 'XYZ input requires --frag-charges'
      return
    end if
    call parse_integer_spec(charge_spec,data%nat,charges,ok,error_message)
    if (.not.ok) return
    if (size(charges) /= data%nfrag) then
      error_message = '--frag-charges must contain exactly one charge per fragment'
      return
    end if
    if (sum(charges) /= data%charge) then
      error_message = 'fragment charges do not sum to the XYZ total charge'
      return
    end if
    call move_alloc(charges,data%fragment_charge)
    data%format='xyz'
  end subroutine read_xyz

  subroutine read_gaussian(filename,data,error_message)
    character(len=*),intent(in) :: filename
    type(eda_input_data),intent(out) :: data
    character(len=:),allocatable,intent(out) :: error_message
    character(len=line_len),allocatable :: lines(:)
    character(len=256) :: tokens(32)
    integer,allocatable :: charge_spin(:),raw_frag(:)
    integer :: nlines,io,idx,route_end,title_end,geom_start,geom_end
    integer :: i,ntok,z,coord_start,frag_value,npairs
    logical :: ok,has_fragment
    real(wp) :: x,y,zcoord

    error_message=''
    call read_lines(filename,lines,nlines,error_message)
    if (len(error_message)>0) return

    idx=1
    do while(idx<=nlines)
      if (starts_with(adjustl(lines(idx)),'#')) exit
      idx=idx+1
    end do
    if (idx>nlines) then
      error_message='Gaussian input has no route section beginning with #'
      return
    end if
    route_end=idx
    do while(route_end<=nlines)
      if (len_trim(lines(route_end))==0) exit
      route_end=route_end+1
    end do
    idx=route_end
    call skip_blank(lines,nlines,idx)
    if (idx>nlines) then
      error_message='Gaussian input is missing the title and molecule specification'
      return
    end if
    title_end=idx
    do while(title_end<=nlines)
      if (len_trim(lines(title_end))==0) exit
      title_end=title_end+1
    end do
    idx=title_end
    call skip_blank(lines,nlines,idx)
    if (idx>nlines) then
      error_message='Gaussian input is missing the charge/multiplicity line'
      return
    end if

    call parse_integer_line(lines(idx),64,charge_spin,ok)
    if (.not.ok .or. size(charge_spin)<2 .or. mod(size(charge_spin),2)/=0) then
      error_message='invalid Gaussian charge/multiplicity line: '//trim(lines(idx))
      return
    end if
    data%charge=charge_spin(1)
    data%spin=charge_spin(2)
    npairs=size(charge_spin)/2
    geom_start=idx+1
    call skip_blank(lines,nlines,geom_start)
    geom_end=geom_start
    do while(geom_end<=nlines)
      if (len_trim(lines(geom_end))==0) exit
      if (lowercase(trim(lines(geom_end)))=='--link1--') exit
      geom_end=geom_end+1
    end do
    data%nat=geom_end-geom_start
    if (data%nat<1) then
      error_message='Gaussian Cartesian geometry is empty'
      return
    end if
    allocate(data%at(data%nat),data%xyz(3,data%nat),raw_frag(data%nat))

    do i=1,data%nat
      call whitespace_tokens(lines(geom_start+i-1),tokens,ntok)
      if (ntok<4) then
        error_message='invalid Gaussian Cartesian line '//integer_string(i)
        return
      end if
      z=gaussian_atom_to_z(tokens(1))
      if (z<1) then
        error_message='unsupported Gaussian atom label on line '//integer_string(i)//': '//trim(tokens(1))
        return
      end if
      call extract_gaussian_fragment(tokens(1),frag_value,has_fragment)
      if (.not.has_fragment .or. frag_value<1) then
        error_message='every Gaussian atom must define Fragment=N; atom '//integer_string(i)//' does not'
        return
      end if

      coord_start=2
      if (ntok>=5 .and. is_integer_token(tokens(2))) then
        if (is_real_token(tokens(3)) .and. is_real_token(tokens(4)) .and. is_real_token(tokens(5))) coord_start=3
      end if
      if (ntok<coord_start+2) then
        error_message='incomplete Gaussian Cartesian coordinates on atom '//integer_string(i)
        return
      end if
      read(tokens(coord_start),*,iostat=io) x
      if (io==0) read(tokens(coord_start+1),*,iostat=io) y
      if (io==0) read(tokens(coord_start+2),*,iostat=io) zcoord
      if (io/=0) then
        error_message='symbolic or invalid Gaussian coordinate on atom '//integer_string(i)
        return
      end if
      data%at(i)=z
      data%xyz(:,i)=[x,y,zcoord]*angstrom_to_bohr
      raw_frag(i)=frag_value
    end do

    allocate(data%fragment(data%nat),source=raw_frag)
    call normalize_fragments(data%fragment,data%nfrag,error_message)
    if (len(error_message)>0) return
    if (data%nfrag<2) then
      error_message='EDA-GFNFF requires at least two Gaussian fragments'
      return
    end if

    if (npairs /= data%nfrag+1) then
      error_message='Gaussian fragment input requires total charge/spin followed by one charge/spin pair per fragment'
      return
    end if
    allocate(data%fragment_charge(data%nfrag))
    do i=1,data%nfrag
      data%fragment_charge(i)=charge_spin(2*i+1)
    end do
    if (sum(data%fragment_charge)/=data%charge) then
      error_message='Gaussian fragment charges do not sum to the total charge'
      return
    end if
    data%format='gaussian'
  end subroutine read_gaussian

  subroutine parse_xyz_comment(comment,charge,spin)
    character(len=*),intent(in) :: comment
    integer,intent(out) :: charge,spin
    logical :: found,ok
    integer :: value

    charge=0
    spin=1
    call extract_key_integer(comment,'charge',value,found)
    if (.not.found) call extract_key_integer(comment,'chrg',value,found)
    if (found) then
      charge=value
    else
      call parse_leading_charge_spin(comment,charge,spin,ok)
    end if
    call extract_key_integer(comment,'spin',value,found)
    if (.not.found) call extract_key_integer(comment,'multiplicity',value,found)
    if (found .and. value>0) spin=value
  end subroutine parse_xyz_comment

  subroutine parse_leading_charge_spin(text,charge,spin,ok)
    character(len=*),intent(in) :: text
    integer,intent(inout) :: charge,spin
    logical,intent(out) :: ok
    character(len=256) :: tokens(32)
    integer :: ntok,io,q,m

    ok=.false.
    call whitespace_tokens(text,tokens,ntok)
    if (ntok<2) return
    read(tokens(1),*,iostat=io) q
    if (io/=0) return
    read(tokens(2),*,iostat=io) m
    if (io/=0 .or. m<1) return
    charge=q; spin=m; ok=.true.
  end subroutine parse_leading_charge_spin

  subroutine extract_key_integer(text,key,value,found)
    character(len=*),intent(in) :: text,key
    integer,intent(out) :: value
    logical,intent(out) :: found
    character(len=:),allocatable :: low,tail
    integer :: p,eq,i,j,io

    found=.false.; value=0
    low=lowercase(text)
    p=index(low,lowercase(key))
    if (p==0) return
    i=p+len_trim(key)
    if (i>len_trim(low)) return
    eq=index(low(i:),'=')
    if (eq==0) return
    eq=i+eq-1
    i=eq+1
    do while(i<=len_trim(text))
      if (text(i:i)/=' ' .and. iachar(text(i:i))/=9) exit
      i=i+1
    end do
    j=i
    if (j<=len_trim(text)) then
      if (text(j:j)=='+' .or. text(j:j)=='-') j=j+1
    end if
    do while(j<=len_trim(text))
      if (text(j:j)<'0' .or. text(j:j)>'9') exit
      j=j+1
    end do
    if (j<=i) return
    tail=text(i:j-1)
    read(tail,*,iostat=io) value
    found=io==0
  end subroutine extract_key_integer

  subroutine parse_fragment_token(token,value,found)
    character(len=*),intent(in) :: token
    integer,intent(out) :: value
    logical,intent(out) :: found
    character(len=:),allocatable :: low,tail
    integer :: p,io

    found=.false.; value=0
    low=lowercase(trim(token))
    p=index(low,'frag=')
    if (p==1) then
      tail=token(6:)
    else if (index(low,'fragment=')==1) then
      tail=token(10:)
    else
      tail=trim(token)
    end if
    read(tail,*,iostat=io) value
    found=io==0
  end subroutine parse_fragment_token

  subroutine extract_gaussian_fragment(atom_label,value,found)
    character(len=*),intent(in) :: atom_label
    integer,intent(out) :: value
    logical,intent(out) :: found
    character(len=:),allocatable :: low,tail
    integer :: p,eq,i,j,io

    found=.false.; value=0
    low=lowercase(atom_label)
    p=index(low,'fragment')
    if (p==0) return
    eq=index(low(p:),'=')
    if (eq==0) return
    eq=p+eq-1
    i=eq+1
    j=i
    if (j<=len_trim(atom_label)) then
      if (atom_label(j:j)=='+' .or. atom_label(j:j)=='-') j=j+1
    end if
    do while(j<=len_trim(atom_label))
      if (atom_label(j:j)<'0' .or. atom_label(j:j)>'9') exit
      j=j+1
    end do
    if (j<=i) return
    tail=atom_label(i:j-1)
    read(tail,*,iostat=io) value
    found=io==0
  end subroutine extract_gaussian_fragment

  integer function gaussian_atom_to_z(label) result(z)
    character(len=*),intent(in) :: label
    character(len=:),allocatable :: head
    integer :: p

    head=trim(label)
    p=index(head,'(')
    if (p>0) head=head(:p-1)
    p=index(head,'-')
    if (p>1) head=head(:p-1)
    z=element_to_z(head)
  end function gaussian_atom_to_z

  integer function element_to_z(symbol) result(z)
    character(len=*),intent(in) :: symbol
    character(len=:),allocatable :: cleaned
    integer :: i,io,numeric

    z=0
    read(symbol,*,iostat=io) numeric
    if (io==0) then
      if (numeric>=1 .and. numeric<=118) z=numeric
      return
    end if
    cleaned=nice_element(symbol)
    select case (lowercase(trim(cleaned)))
    case ('d','t')
      z=1; return
    end select
    do i=1,118
      if (lowercase(trim(pse(i)))==lowercase(trim(cleaned))) then
        z=i; return
      end if
    end do
  end function element_to_z

  function nice_element(text) result(out)
    character(len=*),intent(in) :: text
    character(len=:),allocatable :: out
    character(len=len_trim(text)) :: tmp
    integer :: i,n

    tmp=adjustl(trim(text)); n=len_trim(tmp)
    do i=1,n
      if ((tmp(i:i)>='0'.and.tmp(i:i)<='9') .or. index('*_+',tmp(i:i))>0) tmp(i:i)=' '
    end do
    out=adjustl(trim(tmp))
  end function nice_element

  subroutine normalize_fragments(fragment,nfrag,error_message)
    integer,intent(inout) :: fragment(:)
    integer,intent(out) :: nfrag
    character(len=:),allocatable,intent(out) :: error_message
    integer,allocatable :: unique(:)
    integer :: i,j,k,tmp
    logical :: seen

    error_message=''; nfrag=0
    if (any(fragment<1)) then
      error_message='fragment ids must be positive integers'
      return
    end if
    allocate(unique(size(fragment)))
    do i=1,size(fragment)
      seen=.false.
      do j=1,nfrag
        if (unique(j)==fragment(i)) then
          seen=.true.; exit
        end if
      end do
      if (.not.seen) then
        nfrag=nfrag+1; unique(nfrag)=fragment(i)
      end if
    end do
    do i=1,nfrag-1
      do j=i+1,nfrag
        if (unique(j)<unique(i)) then
          tmp=unique(i); unique(i)=unique(j); unique(j)=tmp
        end if
      end do
    end do
    do i=1,size(fragment)
      do k=1,nfrag
        if (fragment(i)==unique(k)) then
          fragment(i)=k; exit
        end if
      end do
    end do
  end subroutine normalize_fragments

  subroutine parse_integer_spec(spec,max_count,values,ok,error_message)
    character(len=*),intent(in) :: spec
    integer,intent(in) :: max_count
    integer,allocatable,intent(out) :: values(:)
    logical,intent(out) :: ok
    character(len=:),allocatable,intent(out) :: error_message
    integer,allocatable :: buffer(:)
    integer :: unit,io,n
    character(len=line_len) :: line
    logical :: is_file

    allocate(buffer(max_count)); n=0; ok=.true.; error_message=''
    inquire(file=trim(spec),exist=is_file)
    if (is_file) then
      open(newunit=unit,file=trim(spec),status='old',action='read',iostat=io)
      if (io/=0) then
        ok=.false.; error_message='cannot open integer-list file: '//trim(spec); return
      end if
      do
        read(unit,'(a)',iostat=io) line
        if (io/=0) exit
        call append_integer_text(line,buffer,n,max_count,ok)
        if (.not.ok) exit
      end do
      close(unit)
    else
      call append_integer_text(spec,buffer,n,max_count,ok)
    end if
    if (.not.ok .or. n==0) then
      ok=.false.
      if (len(error_message)==0) error_message='invalid or too long integer list: '//trim(spec)
      allocate(values(0)); return
    end if
    allocate(values(n),source=buffer(:n))
  end subroutine parse_integer_spec

  subroutine parse_integer_line(text,max_count,values,ok)
    character(len=*),intent(in) :: text
    integer,intent(in) :: max_count
    integer,allocatable,intent(out) :: values(:)
    logical,intent(out) :: ok
    integer,allocatable :: buffer(:)
    integer :: n

    allocate(buffer(max_count)); n=0; ok=.true.
    call append_integer_text(text,buffer,n,max_count,ok)
    if (.not.ok) then
      allocate(values(0)); return
    end if
    allocate(values(n),source=buffer(:n))
  end subroutine parse_integer_line

  subroutine append_integer_text(text,buffer,n,max_count,ok)
    character(len=*),intent(in) :: text
    integer,intent(inout) :: buffer(:),n
    integer,intent(in) :: max_count
    logical,intent(inout) :: ok
    integer :: i,j,io,value,last
    character(len=128) :: token

    if (.not.ok) return
    last=len_trim(text); i=1
    do while(i<=last)
      do while(i<=last)
        if (.not.is_list_separator(text(i:i))) exit
        i=i+1
      end do
      if (i>last) exit
      j=i
      do while(j<=last)
        if (is_list_separator(text(j:j))) exit
        j=j+1
      end do
      token=''; token=adjustl(text(i:j-1))
      read(token,*,iostat=io) value
      if (io/=0) then
        ok=.false.; return
      end if
      if (n>=max_count) then
        ok=.false.; return
      end if
      n=n+1; buffer(n)=value
      i=j+1
    end do
  end subroutine append_integer_text

  logical function is_list_separator(ch)
    character(len=1),intent(in) :: ch
    is_list_separator = ch==' ' .or. ch==',' .or. ch==';' .or. iachar(ch)==9
  end function is_list_separator

  subroutine whitespace_tokens(line,tokens,ntok)
    character(len=*),intent(in) :: line
    character(len=*),intent(out) :: tokens(:)
    integer,intent(out) :: ntok
    integer :: i,j,last

    tokens=''; ntok=0; last=len_trim(line); i=1
    do while(i<=last)
      do while(i<=last)
        if (line(i:i)/=' ' .and. iachar(line(i:i))/=9) exit
        i=i+1
      end do
      if (i>last) exit
      j=i
      do while(j<=last)
        if (line(j:j)==' ' .or. iachar(line(j:j))==9) exit
        j=j+1
      end do
      if (ntok<size(tokens)) then
        ntok=ntok+1; tokens(ntok)=line(i:j-1)
      end if
      i=j+1
    end do
  end subroutine whitespace_tokens

  subroutine read_lines(filename,lines,nlines,error_message)
    character(len=*),intent(in) :: filename
    character(len=line_len),allocatable,intent(out) :: lines(:)
    integer,intent(out) :: nlines
    character(len=:),allocatable,intent(out) :: error_message
    integer :: unit,io,i
    character(len=line_len) :: line

    error_message=''; nlines=0
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=io)
    if (io/=0) then
      error_message='cannot open Gaussian input: '//trim(filename); return
    end if
    do
      read(unit,'(a)',iostat=io) line
      if (io/=0) exit
      nlines=nlines+1
    end do
    rewind(unit)
    allocate(lines(nlines))
    do i=1,nlines
      read(unit,'(a)') lines(i)
    end do
    close(unit)
  end subroutine read_lines

  subroutine skip_blank(lines,nlines,idx)
    character(len=*),intent(in) :: lines(:)
    integer,intent(in) :: nlines
    integer,intent(inout) :: idx
    do while(idx<=nlines)
      if (len_trim(lines(idx))/=0) exit
      idx=idx+1
    end do
  end subroutine skip_blank

  logical function is_integer_token(token)
    character(len=*),intent(in) :: token
    integer :: value,io
    read(token,*,iostat=io) value
    is_integer_token=io==0
  end function is_integer_token

  logical function is_real_token(token)
    character(len=*),intent(in) :: token
    real(wp) :: value
    integer :: io
    read(token,*,iostat=io) value
    is_real_token=io==0
  end function is_real_token

  logical function starts_with(text,prefix)
    character(len=*),intent(in) :: text,prefix
    integer :: n
    n=len_trim(prefix)
    starts_with=.false.
    if (n==0) then
      starts_with=.true.
    else if (len_trim(text)>=n) then
      starts_with=text(:n)==prefix(:n)
    end if
  end function starts_with

  function lowercase(text) result(out)
    character(len=*),intent(in) :: text
    character(len=len(text)) :: out
    integer :: i,c
    out=text
    do i=1,len(text)
      c=iachar(out(i:i))
      if (c>=iachar('A') .and. c<=iachar('Z')) out(i:i)=achar(c+32)
    end do
  end function lowercase

  function file_extension(filename) result(ext)
    character(len=*),intent(in) :: filename
    character(len=:),allocatable :: ext
    integer :: i,p
    p=0
    do i=len_trim(filename),1,-1
      if (filename(i:i)=='.') then
        p=i; exit
      end if
      if (filename(i:i)=='/' .or. iachar(filename(i:i))==92) exit
    end do
    if (p==0 .or. p==len_trim(filename)) then
      ext=''
    else
      ext=lowercase(filename(p+1:len_trim(filename)))
    end if
  end function file_extension

  function integer_string(value) result(text)
    integer,intent(in) :: value
    character(len=:),allocatable :: text
    character(len=64) :: buffer
    write(buffer,'(i0)') value
    text=trim(buffer)
  end function integer_string

end module eda_gfnff_input
