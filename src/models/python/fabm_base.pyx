cdef extern void c_register_variable(void* pbase, const char* name, const char* units, const char* long_name, int domain, int source, int presence, double initial_value, int* read_index, int* sms_index)
cdef extern void c_unpack_horizontal_cache(void* cache, int* ni, int* ni_hz, int* nread, int* nread_hz, int* nread_scalar, int* nwrite_hz, double** read, double** read_hz, double** read_scalar, double** write_hz)

import importlib

cimport numpy as np
import numpy as np

cdef public void* pnext_base;

domain_interior   = 4
domain_horizontal = 8
domain_scalar     = 16
domain_bottom     = 9
domain_surface    = 10

source_unknown                  =  0
source_do                       =  1
source_do_column                =  2
source_do_horizontal            =  3
source_do_bottom                =  4
source_do_surface               =  5
source_constant                 =  6
source_none                     =  6
source_get_vertical_movement    =  7
source_initialize_state         =  8
source_initialize_surface_state =  9
source_initialize_bottom_state  = 10
source_check_state              = 11
source_check_surface_state      = 12
source_check_bottom_state       = 13
source_get_light_extinction     = 14
source_get_drag                 = 15
source_get_albedo               = 16
source_external                 = 17
source_state                    = 18

presence_internal          = 1
presence_external_required = 2
presence_external_optional = 6

cdef class VariableId:
   cdef str name
   cdef int domain
   cdef int read_index
   cdef int sms_index

cdef class Namespace:
   cdef dict __dict__

cdef class BaseModel:
   cdef void* pbase;
   variables = []
   cdef void* horizontal_cache_source
   cdef dict horizontal_cache

   def __cinit__(self):
      self.pbase = pnext_base
      self.horizontal_cache_source = NULL

   def __dealloc__(self):
      pass

   def register_bottom_state_variable(self, str name, str units, str long_name, *, double initial_value=0.0):
      self.register_variable(name, units, long_name, domain_bottom, source_state, presence_internal, initial_value)

   def register_bottom_dependency(self, str name, str units, str long_name):
      self.register_variable(name, units, long_name, domain_bottom, source_unknown, presence_external_required)

   def register_interior_dependency(self, str name, str units, str long_name):
      self.register_variable(name, units, long_name, domain_interior, source_unknown, presence_external_required)

   def register_variable(self, str name, str units, str long_name, int domain, int source, int presence, double initial_value=0.0):
      cdef VariableId id = VariableId()
      id.domain = domain
      id.name = name
      id.read_index = -1
      id.sms_index = -1
      c_register_variable(self.pbase, name.encode('ascii'), units.encode('ascii'), long_name.encode('ascii'), domain, source, presence, initial_value, &id.read_index, &id.sms_index)
      self.variables.append(id)

   cdef _unpack_horizontal_cache(self, void* cache):
      cdef int ni, ni_hz, nread, nread_hz, nread_scalar, nwrite_hz
      cdef double* read
      cdef double* read_hz
      cdef double* read_scalar
      cdef double* write_hz
      cdef np.ndarray arr_read, arr_read_hz, arr_read_scalar
      cdef VariableId varid

      c_unpack_horizontal_cache(cache, &ni, &ni_hz, &nread, &nread_hz, &nread_scalar, &nwrite_hz, &read, &read_hz, &read_scalar, &write_hz)

      if nread == 0:
         arr_read = None
      elif ni == 1:
         arr_read = np.asarray(<double[:nread:1]> read)
      else:
         arr_read = np.asarray(<double[:nread, :ni:1]> read)

      if nread_hz == 0:
         arr_read_hz = None
         arr_write_hz = None
      elif ni_hz == 1:
         arr_read_hz = np.asarray(<double[:nread_hz:1]> read_hz)
         arr_write_hz = np.asarray(<double[:nwrite_hz:1]> write_hz)
      else:
         arr_read_hz = np.asarray(<double[:nread_hz, :ni_hz:1]> read_hz)
         arr_write_hz = np.asarray(<double[:nwrite_hz, :ni_hz:1]> write_hz)

      if nread_scalar == 0:
         arr_read_scalar = None
      else:
         arr_read_scalar = np.asarray(<double[:nread_scalar:1]> read_scalar)

      self.horizontal_cache = {}
      for varid in self.variables:
         print(varid.name, varid.read_index, varid.sms_index)
         if varid.domain == domain_interior:
            self.horizontal_cache[varid.name] = arr_read[varid.read_index - 1, ...]
         elif varid.domain == domain_bottom or varid.domain == domain_surface:
            self.horizontal_cache[varid.name] = arr_read_hz[varid.read_index - 1, ...]
            if varid.sms_index >= 0:
               self.horizontal_cache[varid.name + '.source'] = arr_write_hz[varid.sms_index - 1, ...]
         elif varid.domain == domain_scalar:
            self.horizontal_cache[varid.name] = arr_read_scalar[varid.read_index - 1, ...]

      self.horizontal_cache_source = cache

cdef public object embedded_python_get_model2(const char* module_name, const char* class_name, void* base):
   pnext = base
   class_name2 = class_name.decode('ascii')
   print(class_name2)
   module_name2 = module_name.decode('ascii')
   print(module_name2)
   mod = importlib.import_module(module_name2)
   cls = getattr(mod, class_name2)
   return cls()

cdef public void embedded_python_do_bottom(object pobject, void* cache):
   cdef BaseModel self = pobject
   if self.horizontal_cache_source != cache:
      self._unpack_horizontal_cache(cache)
   self.do_bottom(self.horizontal_cache)