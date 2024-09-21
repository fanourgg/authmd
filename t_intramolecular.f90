module t_intramolecular_mod

type t_bondtype
   character(len=30) :: type
   double precision, dimension(:), allocatable :: params
   integer                                     :: n !< number of bonds for this bond_type
   integer, dimension(:), allocatable          :: ia, ja
end type t_bondtype

type t_bond
   integer :: n ! total bonds in system
   integer :: ntypes
   type (t_bondtype), dimension(:), allocatable :: type
end type t_bond
!
type t_angletype
   character(len=30) :: type
   double precision, dimension(:), allocatable :: params
   integer                                     :: n !< number of angles for this angle_type
   integer, dimension(:), allocatable          :: ia, ja, ka
end type t_angletype

type t_angle
   integer :: n ! total angles in system
   integer :: ntypes
   type (t_angletype), dimension(:), allocatable :: type
end type t_angle
!
type t_dihedraltype
   character(len=30) :: type
   double precision, dimension(:), allocatable :: params
   integer                                     :: n !< number of dihedrals for this dihedral_type
   integer, dimension(:), allocatable          :: ia, ja, ka, la
end type t_dihedraltype

type t_dihedral
   integer :: n ! total dihedrals in system
   integer :: ntypes
   type (t_dihedraltype), dimension(:), allocatable :: type
end type t_dihedral
!
type t_impropertype
   character(len=30) :: type
   double precision, dimension(:), allocatable :: params
   integer                                     :: n !< number of impropers for this improper_type
   integer, dimension(:), allocatable          :: ia, ja, ka, la
end type t_impropertype

type t_improper
   integer :: n ! total impropers in system
   integer :: ntypes
   type (t_impropertype), dimension(:), allocatable :: type
end type t_improper


type t_intra
   type(t_bond) :: bonds
   type(t_angle) :: angles
   type(t_dihedral) :: dihedrals
   type(t_improper) :: impropers
end type t_intra

type t_intra_envir
   double precision :: ebond, eangle, edihedral, eimproper
   double precision, dimension(3, 3) :: v
end type t_intra_envir


end module t_intramolecular_mod
