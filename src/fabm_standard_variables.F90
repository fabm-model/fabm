module fabm_standard_variables

   ! ====================================================================================================
   ! Data types that contain all metadata needed to describe standard variables.
   ! ====================================================================================================

   type type_bulk_standard_variable
      character(len=256) :: name  = ''
      character(len=64)  :: units = ''
   end type
   
   type type_horizontal_standard_variable
      character(len=256) :: name  = ''
      character(len=64)  :: units = ''
   end type
   
   type type_global_standard_variable
      character(len=256) :: name  = ''
      character(len=64)  :: units = ''
   end type

   ! ====================================================================================================
   ! Standard physical-biogeochemical variables that have a well-defined interpretation and unit.
   ! These variables can be used from biogeochemical models by calling register_dependency with
   ! the name of the desired variable during model initialization.
   ! ====================================================================================================
   ! Names are based on the Standard Name Table from the NetCDF Climate and Forecast (CF) Metadata
   ! Convention. See http://cf-pcmdi.llnl.gov/documents/cf-standard-names/.
   ! In deriving names from the CF convention, the following exceptions are made to account for the fact
   ! that FABM handles both marine and limnic systems and has the water column as default domain:
   ! - "sea_water_" prefix is suppressed
   ! - "_in_sea_water" suffix is suppressed
   ! - instead of the "_at_sea_floor" suffix a "bottom_" prefix is used, analogous to the "surface_"
   !   prefix used in CF.
   ! ====================================================================================================
   
   ! ------------------------------------------------------------------------------------------
   ! Variables defined on the full model domain (e.g., pelagic)
   ! ------------------------------------------------------------------------------------------

   type (type_bulk_standard_variable),parameter :: &
   
      alkalinity_expressed_as_mole_equivalent = type_bulk_standard_variable( &
         'alkalinity_expressed_as_mole_equivalent', &
         'mmol m^-3'), &
      attenuation_coefficient_of_photosynthetic_radiative_flux = type_bulk_standard_variable( &
         'attenuation_coefficient_of_downwelling_photosynthetic_radiative_flux', &
         'm^-1'), &
      cell_thickness = type_bulk_standard_variable( &
         'cell_thickness', &
         'm'), &
      density = type_bulk_standard_variable( &
         'density', &
         'kg m^-3'), &
      downwelling_photosynthetic_radiative_flux = type_bulk_standard_variable( &
         'downwelling_photosynthetic_radiative_flux', &
         'W m^-2'), &
      downwelling_shortwave_flux = type_bulk_standard_variable( &
         'downwelling_shortwave_flux', &
         'W m^-2'), &
      mass_concentration_of_suspended_matter = type_bulk_standard_variable( &
         'mass_concentration_of_suspended_matter', &
         'g m^-3'), &
      mole_concentration_of_ammonium = type_bulk_standard_variable( &
         'mole_concentration_of_ammonium', &
         'mmol m^-3'), &
      mole_concentration_of_dissolved_inorganic_carbon = type_bulk_standard_variable( &
         'mole_concentration_of_dissolved_inorganic_carbon', &
         'mmol m^-3'), &
      mole_concentration_of_dissolved_iron = type_bulk_standard_variable( &
         'mole_concentration_of_dissolved_iron', &
         'umol m^-3'), &
      mole_concentration_of_nitrate = type_bulk_standard_variable( &
         'mole_concentration_of_nitrate', &
         'mmol m^-3'), &
      mole_concentration_of_phosphate = type_bulk_standard_variable( &
         'mole_concentration_of_phosphate', &
         'mmol m^-3'), &
      mole_concentration_of_silicate = type_bulk_standard_variable( &
         'mole_concentration_of_silicate', &
         'mmol m^-3'), &
      ph_reported_on_total_scale = type_bulk_standard_variable( &
         'ph_reported_on_total_scale', &
         '1'), &
      practical_salinity = type_bulk_standard_variable( &
         'practical_salinity', &
         '1e-3'), &
      pressure = type_bulk_standard_variable( &
         'pressure', &
         'dbar'), &
      temperature = type_bulk_standard_variable( &
         'temperature', &
         'degree_Celsius')

   ! ------------------------------------------------------------------------------------------
   ! Variables defined on a horizontal slice of the model domain (e.g., at surface or bottom)
   ! ------------------------------------------------------------------------------------------

   type (type_horizontal_standard_variable),parameter :: &
   
      bottom_depth = type_horizontal_standard_variable( &
         'bottom_depth', &
         'm'), &
      bottom_stress = type_horizontal_standard_variable( &
         'bottom_stress', &
         'Pa'), &
      cloud_area_fraction = type_horizontal_standard_variable( &
         'cloud_area_fraction', &
         '1'), &
      downwelling_photosynthetic_radiative_flux_in_air = type_horizontal_standard_variable( &
         'downwelling_photosynthetic_radiative_flux_in_air', &
         'W m^-2'), &
      downwelling_shortwave_flux_in_air = type_horizontal_standard_variable( &
         'downwelling_shortwave_flux_in_air', &
         'W m^-2'), &
      latitude = type_horizontal_standard_variable( &
         'latitude', &
         'degree_North'), &
      longitude = type_horizontal_standard_variable( &
         'longitude', &
         'degree_East'), &
      mole_fraction_of_carbon_dioxide_in_air = type_horizontal_standard_variable( &
         'mole_fraction_of_carbon_dioxide_in_air', &
         '1e-6'), &
      wind_speed = type_horizontal_standard_variable( &
         'wind_speed', &
         'm s^-1')

   ! ------------------------------------------------------------------------------------------
   ! Global (space-independent) variables
   ! ------------------------------------------------------------------------------------------

   type (type_global_standard_variable),parameter ::  &
   
      number_of_days_since_start_of_the_year = type_global_standard_variable( &
         'number_of_days_since_start_of_the_year', &
         'day')

end module fabm_standard_variables