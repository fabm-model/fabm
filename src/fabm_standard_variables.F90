! ====================================================================================================
! This module declares objects ("identities") for standard physical-biogeochemical variables
! that have a well-defined interpretation and unit.
!
! To use these identities, "use" the fabm_types module, which declares a single "standard_variables"
! variable that has all standard identities as its members.
! (for a list, see the definion of type_standard_variable_collection below)
!
! Biogeochemical models can use these "identity" objects in two ways. First, they can access the value
! of the corresponding variable by registering them as dependency. To do so, call register_dependency
! with the "identity" object as argument.
! Additionally, biogeochemical models can assign a standard identities to their own variables during
! registration. FABM will couple all variables that have been assigned the same identity. Thus, using
! standard identities enables implicit variable coupling.
! ====================================================================================================
! The names of standard variables are based on the Standard Name Table from the NetCDF Climate and
! Forecast (CF) Metadata Convention. See http://cf-pcmdi.llnl.gov/documents/cf-standard-names/.
! In deriving names from the CF convention, the following exceptions are made to account for the fact
! that FABM handles both marine and limnic systems and has the water column as default domain:
! - "sea_water_" prefix is suppressed.
! - "_in_sea_water" suffix is suppressed.
! - instead of the "_at_sea_floor" suffix a "bottom_" prefix is used, analogous to the "surface_"
!   prefix used in CF.
! - the "sea_floor_" prefix is replaced by a "bottom_" prefix.
! ====================================================================================================

module fabm_standard_variables

   private
   
   public type_bulk_standard_variable, type_horizontal_standard_variable, type_global_standard_variable
   public type_standard_variable_collection
   
   ! ====================================================================================================
   ! Data types that contain all metadata needed to describe standard variables.
   ! ====================================================================================================

   type type_bulk_standard_variable
      character(len=256) :: name  = ''    ! Name
      character(len=64)  :: units = ''    ! Units
      character(len=512) :: cf_names = '' ! Comma-separated list of standard names defined in the NetCDF CF convention (http://cf-pcmdi.llnl.gov/documents/cf-standard-names/)
   end type
   
   type type_horizontal_standard_variable
      character(len=256) :: name  = ''    ! Name
      character(len=64)  :: units = ''    ! Units
      character(len=512) :: cf_names = '' ! Comma-separated list of standard names defined in the NetCDF CF convention (http://cf-pcmdi.llnl.gov/documents/cf-standard-names/)
   end type
   
   type type_global_standard_variable
      character(len=256) :: name  = ''    ! Name
      character(len=64)  :: units = ''    ! Units
      character(len=512) :: cf_names = '' ! Comma-separated list of standard names defined in the NetCDF CF convention (http://cf-pcmdi.llnl.gov/documents/cf-standard-names/)
   end type

   ! Single type with all standard variables supported explicitly by FABM.
   ! A single instance of this type is declared in module fabm_types with the name standard_variables.

   type type_standard_variable_collection

      ! Bulk variables
      type (type_bulk_standard_variable) :: &
         alkalinity_expressed_as_mole_equivalent = type_bulk_standard_variable( &
            'alkalinity_expressed_as_mole_equivalent', &
            'mmol m-3', &
            'sea_water_alkalinity_expressed_as_mole_equivalent')
      type (type_bulk_standard_variable) :: &
         attenuation_coefficient_of_photosynthetic_radiative_flux = type_bulk_standard_variable( &
            'attenuation_coefficient_of_photosynthetic_radiative_flux', &
            'm-1', &
            '')
      type (type_bulk_standard_variable) :: &
         cell_thickness = type_bulk_standard_variable( &
            'cell_thickness', &
            'm', &
            'cell_thickness')
      type (type_bulk_standard_variable) :: &
         density = type_bulk_standard_variable( &
            'density', &
            'kg m-3', &
            'sea_water_density')
      type (type_bulk_standard_variable) :: &
         downwelling_photosynthetic_radiative_flux = type_bulk_standard_variable( &
            'downwelling_photosynthetic_radiative_flux', &
            'W m-2', &
            'downwelling_photosynthetic_radiative_flux_in_sea_water')
      type (type_bulk_standard_variable) :: &
         downwelling_shortwave_flux = type_bulk_standard_variable( &
            'downwelling_shortwave_flux', &
            'W m-2', &
            'downwelling_shortwave_flux_in_sea_water')
      type (type_bulk_standard_variable) :: &
         mass_concentration_of_suspended_matter = type_bulk_standard_variable( &
            'mass_concentration_of_suspended_matter', &
            'g m-3', &
            'mass_concentration_of_suspended_matter_in_sea_water,concentration_of_suspended_matter_in_sea_water')
      type (type_bulk_standard_variable) :: &
         mole_concentration_of_ammonium = type_bulk_standard_variable( &
            'mole_concentration_of_ammonium', &
            'mmol m-3', &
            'mole_concentration_of_ammonium_in_sea_water')
      type (type_bulk_standard_variable) :: &
         mole_concentration_of_dissolved_inorganic_carbon = type_bulk_standard_variable( &
            'mole_concentration_of_dissolved_inorganic_carbon', &
            'mmol m-3', &
            'mole_concentration_of_dissolved_inorganic_carbon_in_sea_water')
      type (type_bulk_standard_variable) :: &
         mole_concentration_of_dissolved_iron = type_bulk_standard_variable( &
            'mole_concentration_of_dissolved_iron', &
            'umol m-3', &
            'mole_concentration_of_dissolved_iron_in_sea_water')
      type (type_bulk_standard_variable) :: &
         mole_concentration_of_nitrate = type_bulk_standard_variable( &
            'mole_concentration_of_nitrate', &
            'mmol m-3', &
            'mole_concentration_of_nitrate_in_sea_water')
      type (type_bulk_standard_variable) :: &
         mole_concentration_of_phosphate = type_bulk_standard_variable( &
            'mole_concentration_of_phosphate', &
            'mmol m-3', &
            'mole_concentration_of_phosphate_in_sea_water')
      type (type_bulk_standard_variable) :: &
         mole_concentration_of_silicate = type_bulk_standard_variable( &
            'mole_concentration_of_silicate', &
            'mmol m-3', &
            'mole_concentration_of_silicate_in_sea_water')
      type (type_bulk_standard_variable) :: &
         ph_reported_on_total_scale = type_bulk_standard_variable( &
            'ph_reported_on_total_scale', &
            '1', &
            'sea_water_ph_reported_on_total_scale')
      type (type_bulk_standard_variable) :: &
         practical_salinity = type_bulk_standard_variable( &
            'practical_salinity', &
            '1e-3', &
            'sea_water_practical_salinity')
      type (type_bulk_standard_variable) :: &
         pressure = type_bulk_standard_variable( &
            'pressure', &
            'dbar', &
            'sea_water_pressure')
      type (type_bulk_standard_variable) :: &
         temperature = type_bulk_standard_variable( &
            'temperature', &
            'degree_Celsius', &
            'sea_water_temperature')

      ! Horizontal variables
      type (type_horizontal_standard_variable) :: &
         bottom_depth_below_geoid = type_horizontal_standard_variable( &
            'bottom_depth_below_geoid', &
            'm', &
            'sea_floor_depth_below_geoid,sea_floor_depth')
      type (type_horizontal_standard_variable) :: &
         bottom_stress = type_horizontal_standard_variable( &
            'bottom_stress', &
            'Pa', &
            '')
      type (type_horizontal_standard_variable) :: &
         cloud_area_fraction = type_horizontal_standard_variable( &
            'cloud_area_fraction', &
            '1', &
            'cloud_area_fraction')
      type (type_horizontal_standard_variable) :: &
         surface_downwelling_photosynthetic_radiative_flux = type_horizontal_standard_variable( &
            'surface_downwelling_photosynthetic_radiative_flux', &
            'W m-2', &
            'surface_downwelling_photosynthetic_radiative_flux_in_sea_water')
      type (type_horizontal_standard_variable) :: &
         surface_downwelling_photosynthetic_radiative_flux_in_air = type_horizontal_standard_variable( &
            'surface_downwelling_photosynthetic_radiative_flux_in_air', &
            'W m-2', &
            'surface_downwelling_photosynthetic_radiative_flux_in_air')
      type (type_horizontal_standard_variable) :: &
         surface_downwelling_shortwave_flux = type_horizontal_standard_variable( &
            'surface_downwelling_shortwave_flux', &
            'W m-2', &
            '')
      type (type_horizontal_standard_variable) :: &
         surface_downwelling_shortwave_flux_in_air = type_horizontal_standard_variable( &
            'surface_downwelling_shortwave_flux_in_air', &
            'W m-2', &
            'surface_downwelling_shortwave_flux_in_air,surface_downwelling_shortwave_flux')
      type (type_horizontal_standard_variable) :: &
         latitude = type_horizontal_standard_variable( &
            'latitude', &
            'degree_north', &
            'latitude')
      type (type_horizontal_standard_variable) :: &
         longitude = type_horizontal_standard_variable( &
            'longitude', &
            'degree_east', &
            'longitude')
      type (type_horizontal_standard_variable) :: &
         mole_fraction_of_carbon_dioxide_in_air = type_horizontal_standard_variable( &
            'mole_fraction_of_carbon_dioxide_in_air', &
            '1e-6', &
            'mole_fraction_of_carbon_dioxide_in_air')
      type (type_horizontal_standard_variable) :: &
         wind_speed = type_horizontal_standard_variable( &
            'wind_speed', &
            'm s-1', &
            'wind_speed')

      ! Global variables
      type (type_global_standard_variable) :: &
         number_of_days_since_start_of_the_year = type_global_standard_variable( &
            'number_of_days_since_start_of_the_year', &
            'd', &
            '')

      ! Conserved variables
      type (type_bulk_standard_variable) :: &
         total_carbon = type_bulk_standard_variable( &
            'total_carbon', &
            'mmol m-3', &
            '')
      type (type_bulk_standard_variable) :: &
         total_nitrogen = type_bulk_standard_variable( &
            'total_nitrogen', &
            'mmol m-3', &
            '')
      type (type_bulk_standard_variable) :: &
         total_phosphorus = type_bulk_standard_variable( &
            'total_phosphorus', &
            'mmol m-3', &
            '')
      type (type_bulk_standard_variable) :: &
         total_iron = type_bulk_standard_variable( &
            'total_iron', &
            'umol m-3', &
            '')

   end type type_standard_variable_collection

   ! Collection that contains all standard variables:
   type (type_standard_variable_collection),parameter :: collection = type_standard_variable_collection()

   ! For backward compatibility: individual variables accessible as module-level objects.
   ! Support for these will ultimately disappear; please use the standard_variables
   ! object from fabm_types instead. Then there is no need to access fabm_standard_variables.
   ! For instance: use "standard_variables%temperature" instead of "temperature".

   type (type_bulk_standard_variable), parameter, public :: &
      alkalinity_expressed_as_mole_equivalent = collection%alkalinity_expressed_as_mole_equivalent, &
      attenuation_coefficient_of_photosynthetic_radiative_flux = collection%attenuation_coefficient_of_photosynthetic_radiative_flux, &
      cell_thickness = collection%cell_thickness, &
      density = collection%density, &
      downwelling_photosynthetic_radiative_flux = collection%downwelling_photosynthetic_radiative_flux, &
      downwelling_shortwave_flux = collection%downwelling_shortwave_flux, &
      mass_concentration_of_suspended_matter = collection%mass_concentration_of_suspended_matter, &
      mole_concentration_of_ammonium = collection%mole_concentration_of_ammonium, &
      mole_concentration_of_dissolved_inorganic_carbon = collection%mole_concentration_of_dissolved_inorganic_carbon, &
      mole_concentration_of_dissolved_iron = collection%mole_concentration_of_dissolved_iron, &
      mole_concentration_of_nitrate = collection%mole_concentration_of_nitrate, &
      mole_concentration_of_phosphate = collection%mole_concentration_of_phosphate, &
      mole_concentration_of_silicate = collection%mole_concentration_of_silicate, &
      ph_reported_on_total_scale = collection%ph_reported_on_total_scale, &
      practical_salinity = collection%practical_salinity, &
      pressure = collection%pressure, &
      temperature = collection%temperature

   type (type_horizontal_standard_variable), parameter, public :: &
      bottom_depth = collection%bottom_depth_below_geoid, &
      bottom_depth_below_geoid = collection%bottom_depth_below_geoid, &
      bottom_stress = collection%bottom_stress, &
      cloud_area_fraction = collection%cloud_area_fraction, &
      surface_downwelling_photosynthetic_radiative_flux = collection%surface_downwelling_photosynthetic_radiative_flux, &
      downwelling_photosynthetic_radiative_flux_in_air = collection%surface_downwelling_photosynthetic_radiative_flux_in_air, &
      surface_downwelling_photosynthetic_radiative_flux_in_air = collection%surface_downwelling_photosynthetic_radiative_flux_in_air, &
      surface_downwelling_shortwave_flux = collection%surface_downwelling_shortwave_flux, &
      downwelling_shortwave_flux_in_air = collection%surface_downwelling_shortwave_flux_in_air, &
      surface_downwelling_shortwave_flux_in_air = collection%surface_downwelling_shortwave_flux_in_air, &
      latitude = collection%latitude, &
      longitude = collection%longitude, &
      mole_fraction_of_carbon_dioxide_in_air = collection%mole_fraction_of_carbon_dioxide_in_air, &
      wind_speed = collection%wind_speed

   type (type_global_standard_variable), parameter, public :: &
      number_of_days_since_start_of_the_year = collection%number_of_days_since_start_of_the_year

   type (type_bulk_standard_variable), parameter, public :: &
      total_carbon = collection%total_carbon, &
      total_nitrogen = collection%total_nitrogen, &
      total_phosphorus = collection%total_phosphorus, &
      total_iron = collection%total_iron

   end module fabm_standard_variables