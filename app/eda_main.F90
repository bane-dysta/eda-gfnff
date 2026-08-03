! This file is part of gfnff.
! SPDX-License-Identifier: LGPL-3.0-or-later
!
! Standalone EDA-GFNFF command-line interface.
program eda_gfnff_cli
  use iso_fortran_env, only : wp => real64, stdout => output_unit, stderr => error_unit
#ifdef WITH_OpenMP
  use omp_lib
#endif
  use gfnff_interface, only : gfnff_data
  use eda_gfnff_input, only : eda_input_data, read_eda_input, atomic_symbol
  implicit none

  type(eda_input_data) :: input
  type(gfnff_data) :: calculator
  character(len=4096) :: inputfile,fragment_spec,charge_spec,outputfile,arg,nextarg
  character(len=:),allocatable :: error_message
  real(wp),allocatable :: gradient(:,:)
  real(wp) :: energy
  integer :: io,threads,printlevel,i,nargs,index_arg
  logical :: have_input,have_output

  inputfile=''; fragment_spec=''; charge_spec=''; outputfile=''
  threads=1; printlevel=1; have_input=.false.; have_output=.false.
  nargs=command_argument_count()
  if (nargs==0) call print_help(0)

  index_arg=1
  do while(index_arg<=nargs)
    call get_command_argument(index_arg,arg)
    select case(trim(arg))
    case('-h','--help')
      call print_help(0)
    case('-i','--input')
      call require_next(index_arg,nargs,nextarg,'missing value after --input')
      inputfile=trim(nextarg); have_input=.true.
    case('--frag','--fragments')
      call require_next(index_arg,nargs,nextarg,'missing value after --frag')
      fragment_spec=trim(nextarg)
    case('--frag-charges','--fragment-charges')
      call require_next(index_arg,nargs,nextarg,'missing value after --frag-charges')
      charge_spec=trim(nextarg)
    case('-o','--output')
      call require_next(index_arg,nargs,nextarg,'missing value after --output')
      outputfile=trim(nextarg); have_output=.true.
    case('-T','--threads')
      call require_next(index_arg,nargs,nextarg,'missing value after --threads')
      read(nextarg,*,iostat=io) threads
      if (io/=0 .or. threads<1) call fatal('invalid thread count: '//trim(nextarg))
    case('-v','--verbose')
      printlevel=2
    case('-q','--quiet')
      printlevel=0
    case default
      if (len_trim(arg)>0 .and. arg(1:1)=='-') then
        call fatal('unknown option: '//trim(arg))
      else if (.not.have_input) then
        inputfile=trim(arg); have_input=.true.
      else
        call fatal('unexpected positional argument: '//trim(arg))
      end if
    end select
    index_arg=index_arg+1
  end do

  if (.not.have_input) call fatal('no input file specified')
  call read_eda_input(trim(inputfile),trim(fragment_spec),trim(charge_spec),input,error_message)
  if (len(error_message)>0) call fatal(error_message)
  if (.not.have_output) outputfile=default_extxyz_name(trim(inputfile))
  if (trim(outputfile)==trim(inputfile)) call fatal('extxyz output must differ from the input file')

#ifdef WITH_OpenMP
  call omp_set_num_threads(threads)
#endif

  allocate(calculator%userinput)
  calculator%userinput%fraglist=input%fragment
  allocate(calculator%userinput%fragcharges(input%nfrag))
  calculator%userinput%fragcharges=real(input%fragment_charge,wp)
  calculator%write_topo=.false.
  calculator%do_eda=.true.

  allocate(gradient(3,input%nat),source=0.0_wp)
  call calculator%init(input%nat,input%at,input%xyz,ichrg=input%charge, &
  & printlevel=printlevel,iostat=io)
  if (io/=0) call fatal('GFN-FF topology initialization failed')

  call calculator%singlepoint(input%nat,input%at,input%xyz,energy,gradient, &
  & printlevel=printlevel,iostat=io)
  if (io/=0) call fatal('GFN-FF single-point calculation failed')

  if (.not.allocated(calculator%res%eda%electrostatic)) then
    call fatal('GFN-FF returned no EDA data')
  end if
  if (calculator%res%eda%nfrag/=input%nfrag) then
    call fatal('GFN-FF fragment count differs from the requested partition')
  end if

  call print_header()
  write(stdout,'(1x,a,1x,a)') 'Input:',trim(inputfile)
  write(stdout,'(1x,a,1x,a)') 'Format:',trim(input%format)
  write(stdout,'(1x,a,i0)') 'Atoms: ',input%nat
  write(stdout,'(1x,a,i0)') 'Total charge: ',input%charge
  write(stdout,'(1x,a,i0)') 'Multiplicity: ',input%spin
  write(stdout,'(1x,a,i0)') 'Fragments: ',input%nfrag
  write(stdout,'(/,1x,a)') 'Fragment definition'
  write(stdout,'(1x,a6,2x,a8,2x,a8)') 'Frag','Atoms','Charge'
  do i=1,input%nfrag
    write(stdout,'(1x,i6,2x,i8,2x,i8)') i,count(input%fragment==i),input%fragment_charge(i)
  end do

  write(stdout,'(/,1x,a,1x,es22.14,1x,a)') 'GFN-FF total energy:',energy,'Eh'
  write(stdout,'(1x,a,1x,es22.14,1x,a)') 'GFN-FF gradient norm:',sqrt(sum(gradient**2)),'Eh/a0'

  call print_eda(calculator,input%nfrag)
  call validate_atomic_contributions(calculator,input%nfrag)
  call print_atomic_contributions(calculator,input)
  call write_extxyz(trim(outputfile),calculator,input,error_message)
  if (len(error_message)>0) call fatal(error_message)
  write(stdout,'(/,1x,a,1x,a)') 'Atomic contributions written to:',trim(outputfile)

contains

  subroutine print_eda(calc,nfrag)
    type(gfnff_data),intent(in) :: calc
    integer,intent(in) :: nfrag
    real(wp),parameter :: eh_to_kcal=627.5094740631_wp
    real(wp),parameter :: eh_to_kj=2625.4996394799_wp
    real(wp) :: eel,erep,edisp,ehb,exb,etotal
    real(wp) :: sum_el,sum_rep,sum_disp,sum_hb,sum_xb,sum_three,sum_total
    integer :: a,b

    sum_el=0.0_wp; sum_rep=0.0_wp; sum_disp=0.0_wp
    sum_hb=0.0_wp; sum_xb=0.0_wp
    write(stdout,'(/,1x,a)') 'EDA-GFNFF interfragment decomposition (kcal/mol)'
    write(stdout,'(1x,a5,1x,a5,6(2x,a14))') 'FragA','FragB','Electrostatic','Repulsion', &
    & 'Dispersion','H-bond','X-bond','Total NCI'
    do a=1,nfrag-1
      do b=a+1,nfrag
        eel=calc%res%eda%electrostatic(a,b)
        erep=calc%res%eda%repulsion(a,b)
        edisp=calc%res%eda%dispersion(a,b)
        ehb=calc%res%eda%hydrogen_bond(a,b)
        exb=calc%res%eda%halogen_bond(a,b)
        etotal=eel+erep+edisp+ehb+exb
        sum_el=sum_el+eel; sum_rep=sum_rep+erep; sum_disp=sum_disp+edisp
        sum_hb=sum_hb+ehb; sum_xb=sum_xb+exb
        write(stdout,'(1x,i5,1x,i5,6(2x,f14.6))') a,b,eel*eh_to_kcal,erep*eh_to_kcal, &
        & edisp*eh_to_kcal,ehb*eh_to_kcal,exb*eh_to_kcal,etotal*eh_to_kcal
      end do
    end do
    sum_three=sum_el+sum_rep+sum_disp
    sum_total=sum_three+sum_hb+sum_xb
    write(stdout,'(1x,a)') repeat('-',108)
    write(stdout,'(1x,a11,6(2x,f14.6))') 'All pairs',sum_el*eh_to_kcal,sum_rep*eh_to_kcal, &
    & sum_disp*eh_to_kcal,sum_hb*eh_to_kcal,sum_xb*eh_to_kcal,sum_total*eh_to_kcal
    write(stdout,'(/,1x,a)') 'Total interfragment EDA-GFNFF'
    write(stdout,'(1x,a,1x,es22.14,1x,a)') 'Electrostatic:',sum_el,'Eh'
    write(stdout,'(1x,a,1x,es22.14,1x,a)') 'Repulsion:    ',sum_rep,'Eh'
    write(stdout,'(1x,a,1x,es22.14,1x,a)') 'Dispersion:   ',sum_disp,'Eh'
    write(stdout,'(1x,a,1x,es22.14,1x,a)') 'H-bond:       ',sum_hb,'Eh'
    write(stdout,'(1x,a,1x,es22.14,1x,a)') 'X-bond:       ',sum_xb,'Eh'
    write(stdout,'(1x,a,1x,es22.14,1x,a)') 'Three-term:   ',sum_three,'Eh'
    write(stdout,'(1x,a,1x,es22.14,1x,a)') 'Total NCI:    ',sum_total,'Eh'
    write(stdout,'(1x,a,1x,f18.6,1x,a)') 'Total NCI:    ',sum_total*eh_to_kcal,'kcal/mol'
    write(stdout,'(1x,a,1x,f18.6,1x,a)') 'Total NCI:    ',sum_total*eh_to_kj,'kJ/mol'
    write(stdout,'(/,1x,a)') 'Definition: cross-fragment EEQ pair electrostatics + nonbonded GFN-FF repulsion +'
    write(stdout,'(1x,a)') 'dispersion + GFN-FF hydrogen-bond and halogen-bond corrections.'
  end subroutine print_eda

  subroutine validate_atomic_contributions(calc,nfrag)
    type(gfnff_data),intent(in) :: calc
    integer,intent(in) :: nfrag

    call check_atomic_term('electrostatic',nfrag,calc%res%eda%electrostatic, &
    & calc%res%eda%atom_electrostatic)
    call check_atomic_term('repulsion',nfrag,calc%res%eda%repulsion, &
    & calc%res%eda%atom_repulsion)
    call check_atomic_term('dispersion',nfrag,calc%res%eda%dispersion, &
    & calc%res%eda%atom_dispersion)
    call check_atomic_term('hydrogen-bond',nfrag,calc%res%eda%hydrogen_bond, &
    & calc%res%eda%atom_hydrogen_bond)
    call check_atomic_term('halogen-bond',nfrag,calc%res%eda%halogen_bond, &
    & calc%res%eda%atom_halogen_bond)
  end subroutine validate_atomic_contributions

  subroutine check_atomic_term(name,nfrag,matrix,atoms)
    character(len=*),intent(in) :: name
    integer,intent(in) :: nfrag
    real(wp),intent(in) :: matrix(:,:),atoms(:)
    real(wp),parameter :: tolerance=5.0e-11_wp
    real(wp) :: matrix_sum,atom_sum
    integer :: a,b

    matrix_sum=0.0_wp
    do a=1,nfrag-1
      do b=a+1,nfrag
        matrix_sum=matrix_sum+matrix(a,b)
      end do
    end do
    atom_sum=sum(atoms)
    if (abs(matrix_sum-atom_sum)>tolerance*max(1.0_wp,abs(matrix_sum))) then
      call fatal('atomic '//trim(name)//' contributions do not close to the fragment-pair sum')
    end if
  end subroutine check_atomic_term

  subroutine print_atomic_contributions(calc,data)
    type(gfnff_data),intent(in) :: calc
    type(eda_input_data),intent(in) :: data
    real(wp),parameter :: eh_to_kcal=627.5094740631_wp
    real(wp) :: eel,erep,edisp,ehb,exb,three,total
    real(wp) :: sum_el,sum_rep,sum_disp,sum_hb,sum_xb,sum_three,sum_total
    integer :: i

    write(stdout,'(/,1x,a)') 'Multiwfn-style atomic contributions (kcal/mol)'
    write(stdout,'(1x,a6,1x,a4,1x,a5,7(2x,a13))') 'Atom','Elem','Frag','Electrostatic', &
    & 'Repulsion','Dispersion','H-bond','X-bond','Three-term','Total NCI'
    do i=1,data%nat
      eel=calc%res%eda%atom_electrostatic(i)
      erep=calc%res%eda%atom_repulsion(i)
      edisp=calc%res%eda%atom_dispersion(i)
      ehb=calc%res%eda%atom_hydrogen_bond(i)
      exb=calc%res%eda%atom_halogen_bond(i)
      three=eel+erep+edisp
      total=three+ehb+exb
      write(stdout,'(1x,i6,1x,a4,1x,i5,7(2x,f13.6))') i,trim(atomic_symbol(data%at(i))), &
      & data%fragment(i),eel*eh_to_kcal,erep*eh_to_kcal,edisp*eh_to_kcal, &
      & ehb*eh_to_kcal,exb*eh_to_kcal,three*eh_to_kcal,total*eh_to_kcal
    end do

    sum_el=sum(calc%res%eda%atom_electrostatic)
    sum_rep=sum(calc%res%eda%atom_repulsion)
    sum_disp=sum(calc%res%eda%atom_dispersion)
    sum_hb=sum(calc%res%eda%atom_hydrogen_bond)
    sum_xb=sum(calc%res%eda%atom_halogen_bond)
    sum_three=sum_el+sum_rep+sum_disp
    sum_total=sum_three+sum_hb+sum_xb
    write(stdout,'(1x,a)') repeat('-',126)
    write(stdout,'(1x,a17,7(2x,f13.6))') 'Atomic sum',sum_el*eh_to_kcal, &
    & sum_rep*eh_to_kcal,sum_disp*eh_to_kcal,sum_hb*eh_to_kcal, &
    & sum_xb*eh_to_kcal,sum_three*eh_to_kcal,sum_total*eh_to_kcal
  end subroutine print_atomic_contributions

  subroutine write_extxyz(filename,calc,data,error_message)
    character(len=*),intent(in) :: filename
    type(gfnff_data),intent(in) :: calc
    type(eda_input_data),intent(in) :: data
    character(len=:),allocatable,intent(out) :: error_message
    real(wp),parameter :: bohr_to_angstrom=0.529177249_wp
    real(wp),parameter :: eh_to_kcal=627.5094740631_wp
    real(wp) :: eel,erep,edisp,ehb,exb,three,total,total_nci
    integer :: unit,io,i
    character(len=64) :: total_text

    error_message=''
    open(newunit=unit,file=filename,status='replace',action='write',iostat=io)
    if (io/=0) then
      error_message='cannot create extxyz output: '//trim(filename)
      return
    end if

    total_nci=sum(calc%res%eda%atom_electrostatic)+sum(calc%res%eda%atom_repulsion)+ &
    & sum(calc%res%eda%atom_dispersion)+sum(calc%res%eda%atom_hydrogen_bond)+ &
    & sum(calc%res%eda%atom_halogen_bond)
    write(total_text,'(es24.16)') total_nci*eh_to_kcal
    write(unit,'(i0)') data%nat
    write(unit,'(a)') &
    & 'Properties=species:S:1:pos:R:3:fragment:I:1:eda_electrostatic:R:1:'// &
    & 'eda_repulsion:R:1:eda_dispersion:R:1:eda_hbond:R:1:eda_xbond:R:1:'// &
    & 'eda_three_term:R:1:eda_total_nci:R:1 energy_unit="kcal/mol" '// &
    & 'atomic_contribution="Multiwfn pair-halving convention; HB/XB equally split over centers" '// &
    & 'charge='//trim(integer_text(data%charge))//' multiplicity='//trim(integer_text(data%spin))// &
    & ' nfrag='//trim(integer_text(data%nfrag))//' total_nci='//trim(adjustl(total_text))// &
    & ' source="EDA-GFNFF 0.3"'

    do i=1,data%nat
      eel=calc%res%eda%atom_electrostatic(i)*eh_to_kcal
      erep=calc%res%eda%atom_repulsion(i)*eh_to_kcal
      edisp=calc%res%eda%atom_dispersion(i)*eh_to_kcal
      ehb=calc%res%eda%atom_hydrogen_bond(i)*eh_to_kcal
      exb=calc%res%eda%atom_halogen_bond(i)*eh_to_kcal
      three=eel+erep+edisp
      total=three+ehb+exb
      write(unit,'(a,3(1x,f18.10),1x,i0,7(1x,es20.12))') trim(atomic_symbol(data%at(i))), &
      & data%xyz(:,i)*bohr_to_angstrom,data%fragment(i),eel,erep,edisp,ehb,exb,three,total
    end do
    close(unit)
  end subroutine write_extxyz

  function default_extxyz_name(filename) result(output)
    character(len=*),intent(in) :: filename
    character(len=:),allocatable :: output
    integer :: slash,dot

    slash=max(index(filename,'/',back=.true.),index(filename,achar(92),back=.true.))
    dot=index(filename,'.',back=.true.)
    if (dot>slash+1) then
      output=filename(:dot-1)//'.eda.extxyz'
    else
      output=trim(filename)//'.eda.extxyz'
    end if
  end function default_extxyz_name

  function integer_text(value) result(text)
    integer,intent(in) :: value
    character(len=32) :: text
    write(text,'(i0)') value
  end function integer_text

  subroutine require_next(position,total,value,message)
    integer,intent(inout) :: position
    integer,intent(in) :: total
    character(len=*),intent(out) :: value
    character(len=*),intent(in) :: message
    if (position>=total) call fatal(message)
    position=position+1
    call get_command_argument(position,value)
  end subroutine require_next

  subroutine fatal(message)
    character(len=*),intent(in) :: message
    write(stderr,'(a)') 'eda-gfnff: '//trim(message)
    error stop 2
  end subroutine fatal

  subroutine print_header()
    write(stdout,'(a)') ''
    write(stdout,'(1x,a)') 'EDA-GFNFF 0.3'
    write(stdout,'(1x,a)') repeat('=',60)
  end subroutine print_header

  subroutine print_help(status)
    integer,intent(in) :: status
    write(stdout,'(a)') 'Usage:'
    write(stdout,'(a)') '  eda-gfnff complex.xyz --frag FRAGMENTS --frag-charges CHARGES [options]'
    write(stdout,'(a)') '  eda-gfnff complex.gjf [options]'
    write(stdout,'(a)') '  eda-gfnff complex.com [options]'
    write(stdout,'(a)') ''
    write(stdout,'(a)') 'XYZ input:'
    write(stdout,'(a)') '  Comment line: "0 1", "charge=-1", optionally spin=2/multiplicity=2.'
    write(stdout,'(a)') '  Fragment ids: atom-line column 5, or --frag with a list/file.'
    write(stdout,'(a)') '  Fragment charges: --frag-charges with a list/file in ascending fragment-id order.'
    write(stdout,'(a)') '  Example: eda-gfnff dimer.xyz --frag "1,1,1,2,2" --frag-charges "0,0"'
    write(stdout,'(a)') ''
    write(stdout,'(a)') 'Gaussian input:'
    write(stdout,'(a)') '  Cartesian atoms must use Element(Fragment=N). The charge line must contain'
    write(stdout,'(a)') '  total charge/spin followed by one charge/spin pair per fragment, e.g.'
    write(stdout,'(a)') '  0 1  0 1  0 1'
    write(stdout,'(a)') ''
    write(stdout,'(a)') 'Options:'
    write(stdout,'(a)') '  -i, --input FILE'
    write(stdout,'(a)') '  --frag LIST_OR_FILE'
    write(stdout,'(a)') '  --frag-charges LIST_OR_FILE'
    write(stdout,'(a)') '  -o, --output FILE       extxyz output (default: INPUT.eda.extxyz)'
    write(stdout,'(a)') '  -T, --threads N'
    write(stdout,'(a)') '  -v, --verbose'
    write(stdout,'(a)') '  -q, --quiet'
    write(stdout,'(a)') '  -h, --help'
    if (status==0) stop
    error stop status
  end subroutine print_help

end program eda_gfnff_cli
