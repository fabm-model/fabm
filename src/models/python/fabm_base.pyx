cdef extern void c_register_variable(void* pbase, const char* name, const char* units, const char* long_name, int domain, int source, int presence, double initial_value, int* read_index, int* write_index, int* sms_index) nogil
cdef extern void c_unpack_horizontal_cache(void* cache, int* ni, int* ni_hz, int* nread, int* nread_hz, int* nread_scalar, int* nwrite_hz, double** read, double** read_hz, double** read_scalar, double** write_hz) nogil
cdef extern void c_unpack_interior_cache(void* cache, int* ni, int* ni_hz, int* nread, int* nread_hz, int* nread_scalar, int* nwrite, double** read, double** read_hz, double** read_scalar, double** write) nogil
cdef extern void c_log_message(void* pbase, const char* msg) nogil

import importlib
import traceback

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
   cdef int write_index

cdef class Namespace:
   cdef dict __dict__

cdef class Logger:
   cdef void* pbase

   def __cinit__(self, BaseModel owner):
      self.pbase = owner.pbase

   def write(self, s):
      s_strip = s.rstrip()
      if s_strip:
         c_log_message(self.pbase, s_strip.encode("ascii"))


cdef class BaseModel:
   cdef void* pbase;
   variables = []
   cdef void* interior_cache_source
   cdef void* horizontal_cache_source
   cdef dict interior_cache
   cdef dict horizontal_cache
   cdef Logger logger

   def __cinit__(self):
      self.pbase = pnext_base
      self.interior_cache_source = NULL
      self.horizontal_cache_source = NULL
      self.logger = Logger(self)
      self.interior_cache = {}
      self.horizontal_cache = {}

   def __dealloc__(self):
      pass

   def register_interior_state_variable(self, str name, str units, str long_name, *, double initial_value=0.0):
      self.register_variable(name, units, long_name, domain_interior, source_state, presence_internal, initial_value)

   def register_bottom_state_variable(self, str name, str units, str long_name, *, double initial_value=0.0):
      self.register_variable(name, units, long_name, domain_bottom, source_state, presence_internal, initial_value)

   def register_surface_state_variable(self, str name, str units, str long_name, *, double initial_value=0.0):
      self.register_variable(name, units, long_name, domain_surface, source_state, presence_internal, initial_value)

   def register_interior_dependency(self, str name, str units, str long_name):
      self.register_variable(name, units, long_name, domain_interior, source_unknown, presence_external_required)

   def register_bottom_dependency(self, str name, str units, str long_name):
      self.register_variable(name, units, long_name, domain_bottom, source_unknown, presence_external_required)

   def register_surface_dependency(self, str name, str units, str long_name):
      self.register_variable(name, units, long_name, domain_surface, source_unknown, presence_external_required)

   def register_global_dependency(self, str name, str units, str long_name):
      self.register_variable(name, units, long_name, domain_scalar, source_unknown, presence_external_required)

   def register_interior_diagnostic_variable(self, str name, str units, str long_name):
      self.register_variable(name, units, long_name, domain_interior, source_do, presence_internal)

   def register_bottom_diagnostic_variable(self, str name, str units, str long_name):
      self.register_variable(name, units, long_name, domain_bottom, source_do_bottom, presence_internal)

   def register_surface_diagnostic_variable(self, str name, str units, str long_name):
      self.register_variable(name, units, long_name, domain_surface, source_do_surface, presence_internal)

   def register_variable(self, str name, str units, str long_name, int domain, int source, int presence, double initial_value=0.0):
      cdef VariableId id = VariableId()
      id.domain = domain
      id.name = name
      id.read_index = -1
      id.sms_index = -1
      id.write_index = -1
      c_register_variable(self.pbase, name.encode('ascii'), units.encode('ascii'), long_name.encode('ascii'), domain, source, presence, initial_value, &id.read_index, &id.write_index, &id.sms_index)
      self.variables.append(id)

   cdef _unpack_interior_cache(self, void* cache):
      cdef int ni, ni_hz, nread, nread_hz, nread_scalar, nwrite
      cdef double* read
      cdef double* read_hz
      cdef double* read_scalar
      cdef double* write
      cdef np.ndarray arr_read, arr_read_hz, arr_read_scalar

      c_unpack_interior_cache(cache, &ni, &ni_hz, &nread, &nread_hz, &nread_scalar, &nwrite, &read, &read_hz, &read_scalar, &write)

      arr_read = _unpack_cache_array(nread, ni, read)
      arr_read_hz = _unpack_cache_array(nread_hz, ni_hz, read_hz)
      arr_read_scalar = _unpack_cache_array(nread_scalar, 1, read_scalar)
      arr_write = _unpack_cache_array(nwrite, ni, write, writeable=True)
      arr_write = arr_write[1:, ...]    # first field (Fortran index=0) contains zeros

      self._populate_cache(self.interior_cache, arr_read, arr_read_hz, arr_read_scalar, arr_write, None)
      self.interior_cache_source = cache

   cdef _unpack_horizontal_cache(self, void* cache):
      cdef int ni, ni_hz, nread, nread_hz, nread_scalar, nwrite_hz
      cdef double* read
      cdef double* read_hz
      cdef double* read_scalar
      cdef double* write_hz
      cdef np.ndarray arr_read, arr_read_hz, arr_read_scalar

      c_unpack_horizontal_cache(cache, &ni, &ni_hz, &nread, &nread_hz, &nread_scalar, &nwrite_hz, &read, &read_hz, &read_scalar, &write_hz)

      arr_read = _unpack_cache_array(nread, ni, read)
      arr_read_hz = _unpack_cache_array(nread_hz, ni_hz, read_hz)
      arr_read_scalar = _unpack_cache_array(nread_scalar, 1, read_scalar)
      arr_write_hz = _unpack_cache_array(nwrite_hz, ni_hz, write_hz, writeable=True)
      arr_write_hz = arr_write_hz[1:, ...]    # first field (Fortran index=0) contains zeros

      self._populate_cache(self.horizontal_cache, arr_read, arr_read_hz, arr_read_scalar, None, arr_write_hz)
      self.horizontal_cache_source = cache

   cdef _populate_cache(self, cache, np.ndarray arr_read, np.ndarray arr_read_hz, np.ndarray arr_read_scalar, np.ndarray arr_write, np.ndarray arr_write_hz):
      cdef VariableId varid

      cache.clear()
      for varid in self.variables:
         if varid.domain == domain_interior:
            if varid.read_index >= 0:
               cache[varid.name] = arr_read[varid.read_index - 1, ...]
            if varid.write_index >= 0 and arr_write is not None:
               cache[varid.name] = arr_write[varid.write_index - 1, ...]
            if varid.sms_index >= 0 and arr_write is not None:
               cache[varid.name + '.source'] = arr_write[varid.sms_index - 1, ...]
         elif varid.domain == domain_bottom or varid.domain == domain_surface:
            if varid.read_index >= 0:
               cache[varid.name] = arr_read_hz[varid.read_index - 1, ...]
            if varid.write_index >= 0 and arr_write_hz is not None:
               cache[varid.name] = arr_write_hz[varid.write_index - 1, ...]
            if varid.sms_index >= 0 and arr_write_hz is not None:
               cache[varid.name + '.source'] = arr_write_hz[varid.sms_index - 1, ...]
         elif varid.domain == domain_scalar:
            cache[varid.name] = arr_read_scalar[varid.read_index - 1, ...]

   def do(self, cache):
      pass

   def do_bottom(self, cache):
      pass

   def do_surface(self, cache):
      pass

cdef np.ndarray _unpack_cache_array(int n, int ni, double* array, bint writeable=False):
   cdef np.ndarray arr
   if n == 0:
      return None
   elif ni == 1:
      arr = np.asarray(<double[:n:1]> array)
   else:
      arr = np.asarray(<double[:n, :ni:1]> array)
   if not writeable:
      arr.flags["WRITEABLE"] = False
   return arr

cdef public int embedded_python_do(object pobject, void* cache) noexcept:
   cdef BaseModel self = pobject
   try:
      if self.interior_cache_source != cache:
         self._unpack_interior_cache(cache)
      self.do(self.interior_cache)
   except:
      traceback.print_exc(file=self.logger)
      return -1
   return 0

cdef public int embedded_python_do_bottom(object pobject, void* cache) noexcept:
   cdef BaseModel self = pobject
   try:
      if self.horizontal_cache_source != cache:
         self._unpack_horizontal_cache(cache)
      self.do_bottom(self.horizontal_cache)
   except:
      traceback.print_exc(file=self.logger)
      return -1
   return 0

cdef public int embedded_python_do_surface(object pobject, void* cache) noexcept:
   cdef BaseModel self = pobject

   try:
      if self.horizontal_cache_source != cache:
         self._unpack_horizontal_cache(cache)
      self.do_surface(self.horizontal_cache)
   except:
      traceback.print_exc(file=self.logger)
      return -1
   return 0
