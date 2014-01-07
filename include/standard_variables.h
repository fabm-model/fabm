   ! Single type with all standard variables supported explicitly by FABM.
   ! A single instance of this type is declared below with the name standard_variables.
   ! This same instance is made publicly available in module fabm_types.

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
         depth = type_bulk_standard_variable( &
            'depth', &
            'm', &
            'depth')
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
         fractional_saturation_of_oxygen = type_bulk_standard_variable( &
            'fractional_saturation_of_oxygen', &
            '1', &
            'fractional_saturation_of_oxygen_in_sea_water')
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
         mole_concentration_of_carbonate_expressed_as_carbon = type_bulk_standard_variable( &
            'mole_concentration_of_carbonate_expressed_as_carbon', &
            'mmol m-3', &
            'mole_concentration_of_carbonate_expressed_as_carbon_in_sea_water')
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
         bottom_depth = type_horizontal_standard_variable( &
            'bottom_depth', &
            'm', &
            '')
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
         surface_air_pressure = type_horizontal_standard_variable( &
            'surface_air_pressure', &
            'Pa', &
            'surface_air_pressure')
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
         surface_specific_humidity = type_horizontal_standard_variable( &
            'surface_specific_humidity', &
            '1', &
            'surface_specific_humidity')
      type (type_horizontal_standard_variable) :: &
         surface_temperature = type_horizontal_standard_variable( &
            'surface_temperature', &
            'degree_Celsius', &
            'surface_temperature,surface_temperature_where_land,surface_temperature_where_open_sea,surface_temperature_where_snow')
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
         total_iron = type_bulk_standard_variable( &
            'total_iron', &
            'umol m-3', &
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
         total_silicate = type_bulk_standard_variable( &
            'total_silicate', &
            'mmol m-3', &
            '')

   end type type_standard_variable_collection

   ! Collection that contains all standard variables:
   type (type_standard_variable_collection),parameter :: standard_variables = type_standard_variable_collection()

   ! For backward compatibility: individual variables accessible as module-level objects.
   ! Support for these will ultimately disappear; please use the standard_variables
   ! object from fabm_types instead. Then there is no need to access fabm_standard_variables.
   ! For instance: use "standard_variables%temperature" instead of "temperature".

   type (type_bulk_standard_variable), parameter, public :: &
      alkalinity_expressed_as_mole_equivalent = standard_variables%alkalinity_expressed_as_mole_equivalent, &
      attenuation_coefficient_of_photosynthetic_radiative_flux = standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
      cell_thickness = standard_variables%cell_thickness, &
      density = standard_variables%density, &
      depth = standard_variables%depth, &
      downwelling_photosynthetic_radiative_flux = standard_variables%downwelling_photosynthetic_radiative_flux, &
      downwelling_shortwave_flux = standard_variables%downwelling_shortwave_flux, &
      fractional_saturation_of_oxygen = standard_variables%fractional_saturation_of_oxygen, &
      mass_concentration_of_suspended_matter = standard_variables%mass_concentration_of_suspended_matter, &
      mole_concentration_of_ammonium = standard_variables%mole_concentration_of_ammonium, &
      mole_concentration_of_carbonate_expressed_as_carbon = standard_variables%mole_concentration_of_carbonate_expressed_as_carbon, &
      mole_concentration_of_dissolved_inorganic_carbon = standard_variables%mole_concentration_of_dissolved_inorganic_carbon, &
      mole_concentration_of_dissolved_iron = standard_variables%mole_concentration_of_dissolved_iron, &
      mole_concentration_of_nitrate = standard_variables%mole_concentration_of_nitrate, &
      mole_concentration_of_phosphate = standard_variables%mole_concentration_of_phosphate, &
      mole_concentration_of_silicate = standard_variables%mole_concentration_of_silicate, &
      ph_reported_on_total_scale = standard_variables%ph_reported_on_total_scale, &
      practical_salinity = standard_variables%practical_salinity, &
      pressure = standard_variables%pressure, &
      temperature = standard_variables%temperature

   type (type_horizontal_standard_variable), parameter, public :: &
      bottom_depth = standard_variables%bottom_depth, &
      bottom_depth_below_geoid = standard_variables%bottom_depth_below_geoid, &
      bottom_stress = standard_variables%bottom_stress, &
      cloud_area_fraction = standard_variables%cloud_area_fraction, &
      downwelling_photosynthetic_radiative_flux_in_air = standard_variables%surface_downwelling_photosynthetic_radiative_flux_in_air, &
      downwelling_shortwave_flux_in_air = standard_variables%surface_downwelling_shortwave_flux_in_air, &
      latitude = standard_variables%latitude, &
      longitude = standard_variables%longitude, &
      mole_fraction_of_carbon_dioxide_in_air = standard_variables%mole_fraction_of_carbon_dioxide_in_air, &
      surface_air_pressure = standard_variables%surface_air_pressure, &
      surface_downwelling_photosynthetic_radiative_flux = standard_variables%surface_downwelling_photosynthetic_radiative_flux, &
      surface_downwelling_photosynthetic_radiative_flux_in_air = standard_variables%surface_downwelling_photosynthetic_radiative_flux_in_air, &
      surface_downwelling_shortwave_flux = standard_variables%surface_downwelling_shortwave_flux, &
      surface_downwelling_shortwave_flux_in_air = standard_variables%surface_downwelling_shortwave_flux_in_air, &
      surface_specific_humidity = standard_variables%surface_specific_humidity, &
      surface_temperature = standard_variables%surface_temperature, &
      wind_speed = standard_variables%wind_speed

   type (type_global_standard_variable), parameter, public :: &
      number_of_days_since_start_of_the_year = standard_variables%number_of_days_since_start_of_the_year

   type (type_bulk_standard_variable), parameter, public :: &
      total_carbon = standard_variables%total_carbon, &
      total_iron = standard_variables%total_iron, &
      total_nitrogen = standard_variables%total_nitrogen, &
      total_phosphorus = standard_variables%total_phosphorus, &
      total_silicate = standard_variables%total_silicate
