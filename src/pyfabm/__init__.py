import sys
import os
import ctypes
import re
import logging
import enum
from typing import (
    MutableMapping,
    Optional,
    Tuple,
    Iterable,
    Union,
    Callable,
    Mapping,
    Sequence,
    TypeVar,
    List,
    Dict,
)

try:
    import importlib.metadata

    __version__ = importlib.metadata.version("pyfabm")
except ImportError:
    pass

try:
    import numpy as np
except ImportError:
    print("Unable to import NumPy. Please ensure it is installed.")
    sys.exit(1)
import numpy.typing as npt


LOG_CALLBACK = ctypes.CFUNCTYPE(None, ctypes.c_char_p)

name2lib: MutableMapping[str, ctypes.CDLL] = {}


def _find_library(name: str) -> str:
    # Determine potential names of dynamic library.
    libdir, name = os.path.split(name)
    if os.name == "nt":
        names = (f"{name}.dll", f"lib{name}.dll")
    elif os.name == "posix" and sys.platform == "darwin":
        names = (f"lib{name}.dylib",)
    else:
        names = (f"lib{name}.so",)

    # Find FABM dynamic library.
    basedirs = []
    if libdir != "":
        # Library name already includes directory
        basedirs.append(libdir)
    else:
        # Look first in pyfabm directory, then in Python path.
        basedirs.append(os.path.dirname(os.path.abspath(__file__)))
        basedirs.extend(sys.path)
    for basedir in basedirs:
        for possible_name in names:
            path = os.path.join(basedir, possible_name)
            if os.path.isfile(path):
                return path
    raise Exception(
        f"Unable to locate FABM library {name}."
        f" Looked in {basedirs} for a file named {' or '.join(names)}."
    )


def get_lib(name: str) -> ctypes.CDLL:
    if name in name2lib:
        return name2lib[name]

    if os.path.isfile(name):
        path = name
    else:
        path = _find_library(name)

    # Load FABM library.
    lib = ctypes.CDLL(path)
    lib.dtype = ctypes.c_double
    lib.numpy_dtype = np.dtype(lib.dtype).newbyteorder("=")

    # Driver settings (number of spatial dimensions, depth index)
    lib.get_driver_settings.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.get_driver_settings.restype = ctypes.c_void_p

    ndim_c = ctypes.c_int()
    idepthdim_c = ctypes.c_int()
    mask_type = ctypes.c_int()
    variable_bottom_index = ctypes.c_int()
    lib.get_driver_settings(
        ctypes.byref(ndim_c),
        ctypes.byref(idepthdim_c),
        ctypes.byref(mask_type),
        ctypes.byref(variable_bottom_index),
    )
    ndim_int = ndim_c.value
    ndim_hz = ndim_int if idepthdim_c.value == -1 else ndim_int - 1
    lib.ndim_int = ndim_int
    lib.ndim_hz = ndim_hz
    lib.idepthdim = idepthdim_c.value
    lib.mask_type = mask_type.value
    lib.variable_bottom_index = variable_bottom_index.value != 0

    CONTIGUOUS = "CONTIGUOUS"
    arrtype0D = np.ctypeslib.ndpointer(dtype=lib.dtype, ndim=0, flags=CONTIGUOUS)
    arrtype1D = np.ctypeslib.ndpointer(dtype=lib.dtype, ndim=1, flags=CONTIGUOUS)
    arrtypeInterior = np.ctypeslib.ndpointer(
        dtype=lib.dtype, ndim=ndim_int, flags=CONTIGUOUS
    )
    arrtypeHorizontal = np.ctypeslib.ndpointer(
        dtype=lib.dtype, ndim=ndim_hz, flags=CONTIGUOUS
    )
    arrtypeInteriorExt = np.ctypeslib.ndpointer(
        dtype=lib.dtype, ndim=ndim_int + 1, flags=CONTIGUOUS
    )
    arrtypeHorizontalExt = np.ctypeslib.ndpointer(
        dtype=lib.dtype, ndim=ndim_hz + 1, flags=CONTIGUOUS
    )
    arrtypeInteriorExt2 = np.ctypeslib.ndpointer(
        dtype=lib.dtype, ndim=ndim_int + 2, flags=CONTIGUOUS
    )
    arrtypeHorizontalExt2 = np.ctypeslib.ndpointer(
        dtype=lib.dtype, ndim=ndim_hz + 2, flags=CONTIGUOUS
    )

    # Initialization
    lib.create_model.argtypes = [ctypes.c_char_p] + [ctypes.c_int] * ndim_int
    lib.create_model.restype = ctypes.c_void_p
    if ndim_int > 0:
        lib.set_domain_start.argtypes = [ctypes.c_void_p] + [ctypes.c_int] * ndim_int
        lib.set_domain_start.restype = ctypes.c_void_p
        lib.set_domain_stop.argtypes = [ctypes.c_void_p] + [ctypes.c_int] * ndim_int
        lib.set_domain_stop.restype = ctypes.c_void_p

    # Access to model objects (variables, parameters, dependencies, couplings,
    # model instances)
    lib.get_counts.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.get_counts.restype = None
    lib.get_variable_metadata.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    lib.get_variable_metadata.restype = None
    lib.set_variable_save.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
    ]
    lib.set_variable_save.restype = None
    lib.get_variable.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    lib.get_variable.restype = ctypes.c_void_p
    lib.get_parameter_metadata.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.get_parameter_metadata.restype = None
    lib.get_model_metadata.argtypes = [
        ctypes.c_void_p,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.get_model_metadata.restype = None
    lib.get_coupling.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.get_coupling.restype = None
    lib.get_error_state.argtypes = []
    lib.get_error_state.restype = ctypes.c_int
    lib.get_error.argtypes = [ctypes.c_int, ctypes.c_char_p]
    lib.get_error.restype = None
    lib.reset_error_state.argtypes = []
    lib.reset_error_state.restype = None
    lib.set_log_callback.argtypes = [LOG_CALLBACK]
    lib.set_log_callback.restype = None
    if lib.mask_type == 1:
        lib.set_mask.restype = None
        lib.set_mask.argtypes = [
            ctypes.c_void_p,
            np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=ndim_hz, flags=CONTIGUOUS),
        ]
    elif lib.mask_type == 2:
        lib.set_mask.restype = None
        lib.set_mask.argtypes = [
            ctypes.c_void_p,
            np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=ndim_int, flags=CONTIGUOUS),
            np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=ndim_hz, flags=CONTIGUOUS),
        ]
    if lib.variable_bottom_index:
        lib.set_bottom_index.restype = None
        lib.set_bottom_index.argtypes = [
            ctypes.c_void_p,
            np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=ndim_hz, flags=CONTIGUOUS),
        ]

    # Read access to variable attributes
    lib.variable_get_metadata.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    lib.variable_get_metadata.restype = None
    lib.variable_get_background_value.argtypes = [ctypes.c_void_p]
    lib.variable_get_background_value.restype = lib.dtype
    lib.variable_get_missing_value.argtypes = [ctypes.c_void_p]
    lib.variable_get_missing_value.restype = lib.dtype
    lib.variable_get_long_path.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_char_p,
    ]
    lib.variable_get_long_path.restype = None
    lib.variable_get_suitable_masters.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
    lib.variable_get_suitable_masters.restype = ctypes.c_void_p
    lib.variable_get_output.argtypes = [ctypes.c_void_p]
    lib.variable_get_output.restype = ctypes.c_int
    lib.variable_is_required.argtypes = [ctypes.c_void_p]
    lib.variable_is_required.restype = ctypes.c_int
    lib.variable_get_no_river_dilution.argtypes = [ctypes.c_void_p]
    lib.variable_get_no_river_dilution.restype = ctypes.c_int
    lib.variable_get_no_precipitation_dilution.argtypes = [ctypes.c_void_p]
    lib.variable_get_no_precipitation_dilution.restype = ctypes.c_int
    lib.variable_get_property_type.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    lib.variable_get_property_type.restype = ctypes.c_int
    lib.variable_get_real_property.argtypes = [
        ctypes.c_void_p,
        ctypes.c_char_p,
        lib.dtype,
    ]
    lib.variable_get_real_property.restype = lib.dtype
    lib.variable_get_integer_property.argtypes = [
        ctypes.c_void_p,
        ctypes.c_char_p,
        ctypes.c_int,
    ]
    lib.variable_get_integer_property.restype = ctypes.c_int
    lib.variable_get_logical_property.argtypes = [
        ctypes.c_void_p,
        ctypes.c_char_p,
        ctypes.c_int,
    ]
    lib.variable_get_logical_property.restype = ctypes.c_int
    lib.find_standard_variable.argtypes = [ctypes.c_char_p]
    lib.find_standard_variable.restype = ctypes.c_void_p

    # Read/write/reset access to parameters.
    lib.get_real_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    lib.get_real_parameter.restype = lib.dtype
    lib.get_integer_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    lib.get_integer_parameter.restype = ctypes.c_int
    lib.get_logical_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    lib.get_logical_parameter.restype = ctypes.c_int
    lib.get_string_parameter.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
    ]
    lib.get_string_parameter.restype = None
    lib.reset_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.reset_parameter.restype = None
    lib.set_real_parameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p, lib.dtype]
    lib.set_real_parameter.restype = None
    lib.set_integer_parameter.argtypes = [
        ctypes.c_void_p,
        ctypes.c_char_p,
        ctypes.c_int,
    ]
    lib.set_integer_parameter.restype = None
    lib.set_logical_parameter.argtypes = [
        ctypes.c_void_p,
        ctypes.c_char_p,
        ctypes.c_int,
    ]
    lib.set_logical_parameter.restype = None
    lib.set_string_parameter.argtypes = [
        ctypes.c_void_p,
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    lib.set_string_parameter.restype = None

    # Read access to lists of variables (e.g., suitable coupling targets).
    lib.link_list_count.argtypes = [ctypes.c_void_p]
    lib.link_list_count.restype = ctypes.c_int
    lib.link_list_index.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.link_list_index.restype = ctypes.c_void_p
    lib.link_list_finalize.argtypes = [ctypes.c_void_p]
    lib.link_list_finalize.restype = None

    # Routines for sending pointers to state and dependency data.
    lib.link_interior_state_data.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        arrtypeInterior,
    ]
    lib.link_interior_state_data.restype = None
    lib.link_surface_state_data.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        arrtypeHorizontal,
    ]
    lib.link_surface_state_data.restype = None
    lib.link_bottom_state_data.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
        arrtypeHorizontal,
    ]
    lib.link_bottom_state_data.restype = None
    lib.link_interior_data.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
        arrtypeInterior,
    ]
    lib.link_interior_data.restype = None
    lib.link_horizontal_data.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
        arrtypeHorizontal,
    ]
    lib.link_horizontal_data.restype = None
    lib.link_scalar.argtypes = [ctypes.c_void_p, ctypes.c_void_p, arrtype0D]
    lib.link_scalar.restype = None

    # Read access to diagnostic data.
    lib.get_interior_diagnostic_data.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.get_interior_diagnostic_data.restype = ctypes.POINTER(lib.dtype)
    lib.get_horizontal_diagnostic_data.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.get_horizontal_diagnostic_data.restype = ctypes.POINTER(lib.dtype)
    lib.require_data.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
    lib.require_data.restype = None
    lib.get_standard_variable_data.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.get_standard_variable_data.restype = ctypes.POINTER(lib.dtype)

    lib.start.argtypes = [ctypes.c_void_p]
    lib.start.restype = None

    # Routine for retrieving source-sink terms for the interior domain.
    lib.get_sources.argtypes = [
        ctypes.c_void_p,
        lib.dtype,
        arrtypeInteriorExt,
        arrtypeHorizontalExt,
        arrtypeHorizontalExt,
        ctypes.c_int,
        ctypes.c_int,
        arrtypeInterior,
    ]
    lib.get_sources.restype = None
    lib.get_vertical_movement.argtypes = [ctypes.c_void_p, arrtypeInteriorExt]
    lib.get_vertical_movement.restype = None
    lib.get_conserved_quantities.argtypes = [
        ctypes.c_void_p,
        arrtypeHorizontalExt,
        arrtypeInterior,
    ]
    lib.get_conserved_quantities.restype = None
    lib.check_state.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.check_state.restype = ctypes.c_int

    # Routine for getting git repository version information.
    lib.get_version.argtypes = (ctypes.c_int, ctypes.c_char_p)
    lib.get_version.restype = None

    lib.save_settings.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
    lib.save_settings.restype = ctypes.c_void_p

    if ndim_int == 0:
        lib.integrate.argtypes = [
            ctypes.c_void_p,
            ctypes.c_int,
            ctypes.c_int,
            arrtype1D,
            arrtypeInteriorExt,
            arrtypeInteriorExt2,
            lib.dtype,
            ctypes.c_int,
            ctypes.c_int,
            arrtypeInterior,
        ]
        lib.integrate.restype = None

    lib.set_log_callback(log_callback)

    name2lib[name] = lib
    return lib


logger: Optional[logging.Logger] = None


@LOG_CALLBACK
def log_callback(msg: bytes):
    log(msg.decode("ascii"))


def log(msg: str):
    if logger is not None:
        logger.info(msg)
    else:
        print(msg)


INTERIOR_STATE_VARIABLE = 1
SURFACE_STATE_VARIABLE = 2
BOTTOM_STATE_VARIABLE = 3
INTERIOR_DIAGNOSTIC_VARIABLE = 4
HORIZONTAL_DIAGNOSTIC_VARIABLE = 5
CONSERVED_QUANTITY = 6
INTERIOR_DEPENDENCY = 7
HORIZONTAL_DEPENDENCY = 8
SCALAR_DEPENDENCY = 9


class DataType(enum.IntEnum):
    REAL = 1
    INTEGER = 2
    LOGICAL = 3
    STRING = 4


ATTRIBUTE_LENGTH = 256

DISPLAY_MINIMUM = 0
DISPLAY_NORMAL = 1
DISPLAY_ADVANCED = 2

unicodesuperscript = {
    "0": "⁰",
    "1": "¹",
    "2": "²",
    "3": "³",
    "4": "⁴",
    "5": "⁵",
    "6": "⁶",
    "7": "⁷",
    "8": "⁸",
    "9": "⁹",
    "-": "⁻",
}
unicodesubscript = {
    "0": "₀",
    "1": "₁",
    "2": "₂",
    "3": "₃",
    "4": "₄",
    "5": "₅",
    "6": "₆",
    "7": "₇",
    "8": "₈",
    "9": "₉",
}
supnumber = re.compile(r"(?<=\w)(-?\d+)(?=[ \*+\-/]|$)")
supenumber = re.compile(r"(?<=\d)e(-?\d+)(?=[ \*+\-/]|$)")
oldsupminus = re.compile(r"/(\w+)(?:\*\*|\^)(\d+)(?=[ \*+\-/]|$)")
oldsup = re.compile(r"(?<=\w)(?:\*\*|\^)(-?\d+)(?=[ \*+\-/]|$)")
oldsub = re.compile(r"(?<=\w)_(-?\d+)(?=[ \*+\-/]|$)")


def createPrettyUnit(unit: str) -> str:
    def replace_superscript(m: re.Match) -> str:
        return "".join([unicodesuperscript[n] for n in m.group(1)])

    def replace_subscript(m: re.Match) -> str:
        return "".join([unicodesubscript[n] for n in m.group(1)])

    def reple(m: re.Match) -> str:
        return "×10%s" % "".join([unicodesuperscript[n] for n in m.group(1)])

    def reploldminus(m: re.Match) -> str:
        return " %s⁻%s" % (
            m.group(1),
            "".join([unicodesuperscript[n] for n in m.group(2)]),
        )

    # def replold(m):
    #    return u'%s%s' % (m.group(1),u''.join([unicodesuperscript[n] for n in m.group(2)]))
    unit = oldsup.sub(replace_superscript, unit)
    unit = oldsub.sub(replace_subscript, unit)
    unit = supenumber.sub(reple, unit)
    unit = supnumber.sub(replace_superscript, unit)
    # unit = oldsupminus.sub(reploldminus,unit)
    return unit


class FABMException(Exception):
    pass


def hasError() -> bool:
    for lib in name2lib.values():
        if lib.get_error_state() != 0:
            return True
    return False


def getError() -> Optional[str]:
    for lib in name2lib.values():
        if lib.get_error_state() != 0:
            strmessage = ctypes.create_string_buffer(1024)
            lib.get_error(1024, strmessage)
            return strmessage.value.decode("ascii")


NodeValue = TypeVar("NodeValue")
NodeType = Mapping[str, Union["NodeType", NodeValue]]


def printTree(
    root: NodeType, stringmapper: Callable[[NodeValue], str], indent: str = ""
):
    """Print an indented tree of objects, encoded by dictionaries linking the
    names of children to their subtree, or to their object. Objects are finally
    printed as string obtained by calling the provided stringmapper method."""
    for name, item in root.items():
        if isinstance(item, Mapping):
            log(f"{indent}{name}")
            printTree(item, stringmapper, indent + "   ")
        else:
            log(f"{indent}{name} = {stringmapper(item)}")


class VariableProperties:
    def __init__(self, model: "Model", variable_pointer: ctypes.c_void_p):
        self.model = model
        self._pvariable = variable_pointer

    def __getitem__(self, key: str) -> Union[float, int, bool]:
        typecode = self.model.fabm.variable_get_property_type(
            self._pvariable, key.encode("ascii")
        )
        if typecode == DataType.REAL:
            return self.model.fabm.variable_get_real_property(
                self._pvariable, key.encode("ascii"), -1.0
            )
        elif typecode == DataType.INTEGER:
            return self.model.fabm.variable_get_integer_property(
                self._pvariable, key.encode("ascii"), 0
            )
        elif typecode == DataType.LOGICAL:
            return (
                self.model.fabm.variable_get_logical_property(
                    self._pvariable, key.encode("ascii"), 0
                )
                != 0
            )
        raise KeyError


class Variable(object):
    def __init__(
        self,
        model: "Model",
        name: str,
        units: str,
        long_name: str,
        path: Optional[str] = None,
    ):
        self.model = model
        self.name = name
        self.units = units
        self.units_unicode = None if units is None else createPrettyUnit(units)
        self.long_name = long_name or name
        self.path = path or name

    @property
    def long_path(self) -> str:
        return self.long_name

    @property
    def output_name(self) -> str:
        """Name suitable for output (alphanumeric characters and underscores only)"""
        return re.sub(r"\W", "_", self.name)

    @property
    def options(self) -> Optional[Sequence]:
        """Collection of values that this variable can take.
        `None` if the variable is not limited to any particular value."""
        return None

    def __repr__(self) -> str:
        postfix = f"={self.value}" if hasattr(self, "value") else ""
        return f"<{self.name}{postfix}>"


class VariableFromPointer(Variable):
    def __init__(self, model: "Model", variable_pointer: ctypes.c_void_p):
        strname = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strunits = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        model.fabm.variable_get_metadata(
            variable_pointer, ATTRIBUTE_LENGTH, strname, strunits, strlong_name
        )
        name = strname.value.decode("ascii")
        units = strunits.value.decode("ascii")
        long_name = strlong_name.value.decode("ascii")

        super().__init__(model, name, units, long_name)

        self._pvariable = variable_pointer
        self.properties = VariableProperties(self.model, self._pvariable)

    @property
    def long_path(self) -> str:
        """Long model instance name, followed by a slash, followed by long variable name."""
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        self.model.fabm.variable_get_long_path(
            self._pvariable, ATTRIBUTE_LENGTH, strlong_name
        )
        return strlong_name.value.decode("ascii")

    @property
    def missing_value(self) -> float:
        """Value that indicates missing data, for instance, on land."""
        return self.model.fabm.variable_get_missing_value(self._pvariable)

    def getRealProperty(self, name, default=-1.0) -> float:
        return self.model.fabm.variable_get_real_property(
            self._pvariable, name.encode("ascii"), default
        )


class Dependency(VariableFromPointer):
    def __init__(
        self,
        model: "Model",
        variable_pointer: ctypes.c_void_p,
        shape: Tuple[int],
        link_function: Callable[[ctypes.c_void_p, ctypes.c_void_p, np.ndarray], None],
    ):
        super().__init__(model, variable_pointer)
        self._is_set = False
        self._link_function = link_function
        self._shape = shape

    @property
    def value(self) -> Optional[np.ndarray]:
        return None if not self._is_set else self._data

    @value.setter
    def value(self, value: npt.ArrayLike):
        if not self._is_set:
            self.link(np.empty(self._shape, dtype=self.model.fabm.numpy_dtype))
        self._data[...] = value

    def link(self, data: np.ndarray):
        assert data.shape == self._shape, (
            f"{self.name}: shape of provided array {data.shape}"
            f" does not match the shape required {self._shape}"
        )
        self._data = data
        self._link_function(self.model.pmodel, self._pvariable, self._data)
        self._is_set = True

    @property
    def required(self) -> bool:
        return self.model.fabm.variable_is_required(self._pvariable) != 0


class StateVariable(VariableFromPointer):
    def __init__(
        self, model: "Model", variable_pointer: ctypes.c_void_p, data: np.ndarray
    ):
        super().__init__(model, variable_pointer)
        self._data = data

    @property
    def value(self) -> np.ndarray:
        return self._data

    @value.setter
    def value(self, value: npt.ArrayLike):
        self._data[...] = value

    @property
    def background_value(self) -> float:
        return self.model.fabm.variable_get_background_value(self._pvariable)

    @property
    def output(self) -> bool:
        return self.model.fabm.variable_get_output(self._pvariable) != 0

    @property
    def no_river_dilution(self) -> bool:
        return self.model.fabm.variable_get_no_river_dilution(self._pvariable) != 0

    @property
    def no_precipitation_dilution(self) -> bool:
        return (
            self.model.fabm.variable_get_no_precipitation_dilution(self._pvariable) != 0
        )


class DiagnosticVariable(VariableFromPointer):
    def __init__(
        self,
        model: "Model",
        variable_pointer: ctypes.c_void_p,
        index: int,
        horizontal: bool,
    ):
        super().__init__(model, variable_pointer)
        self._data = None
        self._horizontal = horizontal
        self._index = index + 1

    @property
    def value(self) -> Optional[np.ndarray]:
        return self._data

    @property
    def output(self) -> bool:
        """Whether this diagnostic is meant to be included in output by default"""
        return self.model.fabm.variable_get_output(self._pvariable) != 0

    def _set_save(self, value: bool):
        vartype = (
            HORIZONTAL_DIAGNOSTIC_VARIABLE
            if self._horizontal
            else INTERIOR_DIAGNOSTIC_VARIABLE
        )
        self.model.fabm.set_variable_save(
            self.model.pmodel, vartype, self._index, 1 if value else 0
        )

    #: Whether the value of this diagnostic must be calculated, for instance, for output
    save: bool = property(fset=_set_save)


class Parameter(Variable):
    def __init__(
        self,
        model: "Model",
        name: str,
        index: int,
        units: Optional[str] = None,
        long_name: Optional[str] = None,
        type: Optional[DataType] = None,
        has_default: bool = False,
    ):
        super().__init__(model, name, units, long_name)
        self._type = type
        self._index = index + 1
        self._has_default = has_default

    def _get_value(self, *, default: bool = False):
        default = 1 if default else 0
        if self._type == DataType.REAL:
            return self.model.fabm.get_real_parameter(
                self.model.pmodel, self._index, default
            )
        elif self._type == DataType.INTEGER:
            return self.model.fabm.get_integer_parameter(
                self.model.pmodel, self._index, default
            )
        elif self._type == DataType.LOGICAL:
            return (
                self.model.fabm.get_logical_parameter(
                    self.model.pmodel, self._index, default
                )
                != 0
            )
        elif self._type == DataType.STRING:
            result = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
            self.model.fabm.get_string_parameter(
                self.model.pmodel, self._index, default, ATTRIBUTE_LENGTH, result
            )
            return result.value.decode("ascii")

    @property
    def value(self) -> Union[float, int, bool, str]:
        return self._get_value()

    @value.setter
    def value(self, value: Union[float, int, bool, str]):
        settings = self.model._save_state()

        if self._type == DataType.REAL:
            self.model.fabm.set_real_parameter(
                self.model.pmodel, self.name.encode("ascii"), value
            )
        elif self._type == DataType.INTEGER:
            self.model.fabm.set_integer_parameter(
                self.model.pmodel, self.name.encode("ascii"), value
            )
        elif self._type == DataType.LOGICAL:
            self.model.fabm.set_logical_parameter(
                self.model.pmodel, self.name.encode("ascii"), value
            )
        elif self._type == DataType.STRING:
            self.model.fabm.set_string_parameter(
                self.model.pmodel, self.name.encode("ascii"), value.encode("ascii")
            )

        # Update the model configuration
        # (arrays with variables and parameters have changed)
        self.model._update_configuration(settings)

    @property
    def default(self) -> Union[float, int, bool, str, None]:
        """Default value for this parameter (`None` if no default is set)"""
        if not self._has_default:
            return None
        return self._get_value(default=True)

    def reset(self):
        """Reset this parameter to its default value"""
        settings = self.model._save_state()
        self.model.fabm.reset_parameter(self.model.pmodel, self._index)
        self.model._update_configuration(settings)


class StandardVariable:
    def __init__(self, model: "Model", pointer: ctypes.c_void_p):
        self.model = model
        self._pvariable = pointer

    @property
    def value(self) -> np.ndarray:
        horizontal = ctypes.c_int()
        pdata = self.model.fabm.get_standard_variable_data(
            self.model.pmodel, self._pvariable, horizontal
        )
        if horizontal.value == 0:
            shape = self.model.interior_domain_shape
        else:
            shape = self.model.horizontal_domain.shape
        arr = np.ctypeslib.as_array(pdata, shape)
        return arr.view(dtype=self.model.fabm.numpy_dtype)


T = TypeVar("T")


class NamedObjectList(Sequence[T]):
    def __init__(self, *data: Iterable[T]):
        self._data: List[T] = []
        for d in data:
            self._data.extend(d)
        self._lookup: Optional[Dict[str, T]] = None
        self._lookup_ci: Optional[Dict[str, T]] = None

    def __len__(self) -> int:
        return len(self._data)

    def __getitem__(self, key: Union[int, str]) -> T:
        if isinstance(key, str):
            return self.find(key)
        return self._data[key]

    def __contains__(self, key: Union[T, str]) -> bool:
        if isinstance(key, str):
            try:
                self.find(key)
                return True
            except KeyError:
                return False
        return key in self._data

    def index(self, key: Union[T, str], *args) -> int:
        if isinstance(key, str):
            try:
                key = self.find(key)
            except KeyError:
                raise ValueError from None
        return self._data.index(key, *args)

    def __repr__(self) -> str:
        return repr(self._data)

    def __add__(self, other: "NamedObjectList[T]") -> "NamedObjectList[T]":
        return NamedObjectList(self._data, other._data)

    def find(self, name: str, case_insensitive: bool = False) -> T:
        if case_insensitive:
            if self._lookup_ci is None:
                self._lookup_ci = {obj.name.lower(): obj for obj in self._data}
            return self._lookup_ci[name.lower()]
        else:
            if self._lookup is None:
                self._lookup = {obj.name: obj for obj in self._data}
            return self._lookup[name]

    def clear(self):
        self._data.clear()
        self._lookup = None
        self._lookup_ci = None


class Coupling(VariableFromPointer):
    def __init__(self, model: "Model", index: int):
        self._ptarget = ctypes.c_void_p()
        self._psource = ctypes.c_void_p()
        model.fabm.get_coupling(
            model.pmodel,
            index + 1,
            ctypes.byref(self._psource),
            ctypes.byref(self._ptarget),
        )
        super().__init__(model, self._psource)
        self._options = None

    @property
    def value(self) -> Optional[str]:
        if self._psource.value == self._ptarget.value:
            return
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        self.model.fabm.variable_get_long_path(
            self._ptarget, ATTRIBUTE_LENGTH, strlong_name
        )
        return strlong_name.value.decode("ascii")

    @value.setter
    def value(self, value: str):
        log(f"New coupling specified: {value}")
        pass

    @property
    def options(self) -> Sequence[str]:
        if self._options is None:
            self._options: List[str] = []
            plist = self.model.fabm.variable_get_suitable_masters(
                self.model.pmodel, self._psource
            )
            strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
            for i in range(self.model.fabm.link_list_count(plist)):
                variable = self.model.fabm.link_list_index(plist, i + 1)
                self.model.fabm.variable_get_long_path(
                    variable, ATTRIBUTE_LENGTH, strlong_name
                )
                self._options.append(strlong_name.value.decode("ascii"))
            self.model.fabm.link_list_finalize(plist)
        return self._options


class SubModel(object):
    def __init__(self, model: "Model", name: str):
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        iuser = ctypes.c_int()
        model.fabm.get_model_metadata(
            model.pmodel, name.encode("ascii"), ATTRIBUTE_LENGTH, strlong_name, iuser
        )
        self.long_name = strlong_name.value.decode("ascii")
        self.user_created = iuser.value != 0


class Model(object):
    def __init__(
        self,
        path: Union[str, dict] = "fabm.yaml",
        shape: Tuple[int] = (),
        libname: Optional[str] = None,
        start: Optional[Tuple[int]] = None,
        stop: Optional[Tuple[int]] = None,
    ):
        delete = False
        if isinstance(path, dict):
            import tempfile
            import yaml
            import io

            with tempfile.NamedTemporaryFile(
                suffix=".yaml", prefix="fabm", delete=False
            ) as f:
                with io.TextIOWrapper(f, encoding="ascii") as wrapper:
                    yaml.safe_dump(path, wrapper)
                path = f.name
            delete = True

        if libname is None:
            # Pick one of the built-in FABM libraries (0D or 1D)
            ndim = len(shape)
            if ndim > 1:
                raise FABMException(
                    f"Invalid domain shape {shape}."
                    " Domain must have 0 or 1 dimensions."
                )
            libname = {0: "fabm_0d", 1: "fabm_1d"}[ndim]

        self.fabm = get_lib(libname)

        if len(shape) != self.fabm.ndim_int:
            raise FABMException(
                f"Invalid domain shape {shape}."
                f" Domain must have {self.fabm.ndim_int} dimensions."
                f" FABM library: {self.fabm._name}"
            )

        self.fabm.reset_error_state()
        self._cell_thickness = None
        self.pmodel = self.fabm.create_model(path.encode("ascii"), *shape[::-1])
        if hasError():
            raise FABMException(
                f"An error occurred while parsing {path}:\n{getError()}"
            )
        if start is not None:
            self.fabm.set_domain_start(self.pmodel, *[s + 1 for s in start[::-1]])
        if stop is not None:
            self.fabm.set_domain_stop(self.pmodel, *stop[::-1])
        self.interior_domain_shape = tuple(shape)
        self.horizontal_domain_shape = tuple(
            [l for i, l in enumerate(shape) if i != self.fabm.idepthdim]
        )
        if delete:
            os.remove(path)

        # fmt: off
        self.interior_state_variables: NamedObjectList[StateVariable] = NamedObjectList()
        self.surface_state_variables: NamedObjectList[StateVariable] = NamedObjectList()
        self.bottom_state_variables: NamedObjectList[StateVariable] = NamedObjectList()
        self.interior_diagnostic_variables: NamedObjectList[DiagnosticVariable] = NamedObjectList()
        self.horizontal_diagnostic_variables: NamedObjectList[DiagnosticVariable] = NamedObjectList()
        self.conserved_quantities: NamedObjectList[Variable] = NamedObjectList()
        self.parameters: NamedObjectList[Parameter] = NamedObjectList()
        self.interior_dependencies: NamedObjectList[Dependency] = NamedObjectList()
        self.horizontal_dependencies: NamedObjectList[Dependency] = NamedObjectList()
        self.scalar_dependencies: NamedObjectList[Dependency] = NamedObjectList()
        # fmt: on

        self._update_configuration()
        self._mask = None
        self._bottom_index = None

    def link_mask(self, *masks: np.ndarray):
        if self.fabm.mask_type == 0:
            raise FABMException(
                "the underlying FABM library has been compiled without support for masks"
            )
        if len(masks) != self.fabm.mask_type:
            raise FABMException(
                f"link_mask must be provided with {self.fabm.mask_type} masks"
            )
        if len(masks) > 1:
            assert (
                masks[0].shape == self.interior_domain_shape
                and masks[0].dtype == np.intc
                and masks[0].flags["C_CONTIGUOUS"]
            )
        assert (
            masks[-1].shape == self.horizontal_domain_shape
            and masks[-1].dtype == np.intc
            and masks[-1].flags["C_CONTIGUOUS"]
        )
        self._mask = masks
        self.fabm.set_mask(self.pmodel, *self._mask)

    @property
    def mask(self) -> Union[np.ndarray, Sequence[np.ndarray], None]:
        mask = self._mask
        if mask is not None and len(mask) == 1:
            mask = mask[0]
        return mask

    @mask.setter
    def mask(self, values: Union[npt.ArrayLike, Sequence[npt.ArrayLike]]):
        if self.fabm.mask_type == 1:
            values = (values,)
        if len(values) != self.fabm.mask_type:
            raise FABMException(f"mask must be set to {self.fabm.mask_type} values")
        if self._mask is None:
            masks = (np.ones(self.horizontal_domain_shape, dtype=np.intc),)
            if self.fabm.mask_type > 1:
                masks = (np.ones(self.interior_domain_shape, dtype=np.intc),) + masks
            self.link_mask(*masks)
        for value, mask in zip(values, self._mask):
            if value is not mask:
                mask[...] = value

    def link_bottom_index(self, indices: np.ndarray):
        if not self.fabm.variable_bottom_index:
            raise FABMException(
                "the underlying FABM library has been compiled without support for variable bottom indices"
            )
        assert (
            indices.shape == self.horizontal_domain_shape
            and indices.dtype == np.intc
            and indices.flags["C_CONTIGUOUS"]
        )
        self._bottom_index = indices
        self.fabm.set_bottom_index(self.pmodel, self._bottom_index)

    @property
    def bottom_index(self) -> Optional[np.ndarray]:
        return self._bottom_index

    @bottom_index.setter
    def bottom_index(self, indices: npt.ArrayLike):
        if self._bottom_index is None:
            self.link_bottom_index(np.ones(self.horizontal_domain_shape, dtype=np.intc))
        if indices is not self._bottom_index:
            self._bottom_index[...] = indices

    @property
    def state(self) -> np.ndarray:
        if self._state is None:
            raise Exception(
                "State is not available as one contiguous array because interior and"
                " surface/bottom state variables have different shapes."
                " Use interior_state, surface_state, bottom_state attributes instead."
            )
        return self._state

    @state.setter
    def state(self, value: npt.ArrayLike):
        if self._state is None:
            raise Exception(
                "State is not available as one contiguous array because interior and"
                " surface/bottom state variables have different shapes."
                " Use interior_state, surface_state, bottom_state attributes instead."
            )
        if value is not self._state:
            self._state[...] = value

    @property
    def interior_state(self) -> np.ndarray:
        return self._interior_state

    @interior_state.setter
    def interior_state(self, value: npt.ArrayLike):
        if value is not self._interior_state:
            self._interior_state[...] = value

    @property
    def surface_state(self) -> np.ndarray:
        return self._surface_state

    @surface_state.setter
    def surface_state(self, value: npt.ArrayLike):
        if value is not self._surface_state:
            self._surface_state[...] = value

    @property
    def bottom_state(self) -> np.ndarray:
        return self._bottom_state

    @bottom_state.setter
    def bottom_state(self, value: npt.ArrayLike):
        if value is not self._bottom_state:
            self._bottom_state[...] = value

    def link_cell_thickness(self, data: np.ndarray):
        assert (
            data.shape == self.interior_domain_shape
            and data.dtype == self.fabm.numpy_dtype
            and data.flags["C_CONTIGUOUS"]
        )
        self._cell_thickness = data

    def setCellThickness(self, value: npt.ArrayLike):
        if self._cell_thickness is None:
            self.link_cell_thickness(np.empty(self.interior_domain_shape))
        self._cell_thickness[...] = value

    cell_thickness = property(fset=setCellThickness)

    def getSubModel(self, name: str) -> SubModel:
        return SubModel(self, name)

    def save_settings(self, path: str, display: int = DISPLAY_NORMAL):
        """Write model configuration to yaml file"""
        self.fabm.save_settings(self.pmodel, path.encode("ascii"), display)

    def _save_state(self) -> Tuple:
        environment = {}
        for dependency in self.dependencies:
            if dependency.value is not None:
                environment[dependency.name] = dependency.value
        state = {variable.name: variable.value for variable in self.state_variables}
        return environment, state

    def _restore_state(self, data: Tuple):
        environment, state = data
        for dependency in self.dependencies:
            if dependency.name in environment:
                dependency.value = environment[dependency.name]
        for variable in self.state_variables:
            if variable.name in state:
                variable.value = state[variable.name]

    def _update_configuration(self, settings: Optional[Tuple] = None):
        # Get number of model variables per category
        nstate_interior = ctypes.c_int()
        nstate_surface = ctypes.c_int()
        nstate_bottom = ctypes.c_int()
        ndiag_interior = ctypes.c_int()
        ndiag_horizontal = ctypes.c_int()
        ndependencies_interior = ctypes.c_int()
        ndependencies_horizontal = ctypes.c_int()
        ndependencies_scalar = ctypes.c_int()
        nconserved = ctypes.c_int()
        nparameters = ctypes.c_int()
        ncouplings = ctypes.c_int()
        self.fabm.get_counts(
            self.pmodel,
            ctypes.byref(nstate_interior),
            ctypes.byref(nstate_surface),
            ctypes.byref(nstate_bottom),
            ctypes.byref(ndiag_interior),
            ctypes.byref(ndiag_horizontal),
            ctypes.byref(ndependencies_interior),
            ctypes.byref(ndependencies_horizontal),
            ctypes.byref(ndependencies_scalar),
            ctypes.byref(nconserved),
            ctypes.byref(nparameters),
            ctypes.byref(ncouplings),
        )

        # Allocate memory for state variable values, and send ctypes.pointer to
        # this memory to FABM.
        if self.fabm.idepthdim == -1:
            # No depth dimension, so interior and surface/bottom variables have
            # the same shape. Store values for all together in one contiguous array
            self._state = np.empty(
                (nstate_interior.value + nstate_surface.value + nstate_bottom.value,)
                + self.interior_domain_shape,
                dtype=self.fabm.numpy_dtype,
            )
            self._interior_state = self._state[: nstate_interior.value, ...]
            self._surface_state = self._state[
                nstate_interior.value : nstate_interior.value + nstate_surface.value,
                ...,
            ]
            self._bottom_state = self._state[
                nstate_interior.value + nstate_surface.value :, ...
            ]
        else:
            # Surface/bottom variables have one dimension less than interior variables
            # Store values for each variable type in a separate array.
            self._state = None
            self._interior_state = np.empty(
                (nstate_interior.value,) + self.interior_domain_shape,
                dtype=self.fabm.numpy_dtype,
            )
            self._surface_state = np.empty(
                (nstate_surface.value,) + self.horizontal_domain_shape,
                dtype=self.fabm.numpy_dtype,
            )
            self._bottom_state = np.empty(
                (nstate_bottom.value,) + self.horizontal_domain_shape,
                dtype=self.fabm.numpy_dtype,
            )

        # Retrieve variable metadata
        strname = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strunits = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strpath = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        typecode = ctypes.c_int()
        has_default = ctypes.c_int()
        self.interior_state_variables.clear()
        self.surface_state_variables.clear()
        self.bottom_state_variables.clear()
        self.interior_diagnostic_variables.clear()
        self.horizontal_diagnostic_variables.clear()
        self.conserved_quantities.clear()
        self.parameters.clear()
        self.interior_dependencies.clear()
        self.horizontal_dependencies.clear()
        self.scalar_dependencies.clear()
        for i in range(nstate_interior.value):
            values = self._interior_state[i, ...]
            ptr = self.fabm.get_variable(self.pmodel, INTERIOR_STATE_VARIABLE, i + 1)
            self.interior_state_variables._data.append(StateVariable(self, ptr, values))
            self.fabm.link_interior_state_data(self.pmodel, i + 1, values)
        for i in range(nstate_surface.value):
            values = self._surface_state[i, ...]
            ptr = self.fabm.get_variable(self.pmodel, SURFACE_STATE_VARIABLE, i + 1)
            self.surface_state_variables._data.append(StateVariable(self, ptr, values))
            self.fabm.link_surface_state_data(self.pmodel, i + 1, values)
        for i in range(nstate_bottom.value):
            values = self._bottom_state[i, ...]
            ptr = self.fabm.get_variable(self.pmodel, BOTTOM_STATE_VARIABLE, i + 1)
            self.bottom_state_variables._data.append(StateVariable(self, ptr, values))
            self.fabm.link_bottom_state_data(self.pmodel, i + 1, values)
        for i in range(ndiag_interior.value):
            ptr = self.fabm.get_variable(
                self.pmodel, INTERIOR_DIAGNOSTIC_VARIABLE, i + 1
            )
            self.interior_diagnostic_variables._data.append(
                DiagnosticVariable(self, ptr, i, False)
            )
        for i in range(ndiag_horizontal.value):
            ptr = self.fabm.get_variable(
                self.pmodel, HORIZONTAL_DIAGNOSTIC_VARIABLE, i + 1
            )
            self.horizontal_diagnostic_variables._data.append(
                DiagnosticVariable(self, ptr, i, True)
            )
        for i in range(ndependencies_interior.value):
            ptr = self.fabm.get_variable(self.pmodel, INTERIOR_DEPENDENCY, i + 1)
            self.interior_dependencies._data.append(
                Dependency(
                    self, ptr, self.interior_domain_shape, self.fabm.link_interior_data
                )
            )
        for i in range(ndependencies_horizontal.value):
            ptr = self.fabm.get_variable(self.pmodel, HORIZONTAL_DEPENDENCY, i + 1)
            self.horizontal_dependencies._data.append(
                Dependency(
                    self,
                    ptr,
                    self.horizontal_domain_shape,
                    self.fabm.link_horizontal_data,
                )
            )
        for i in range(ndependencies_scalar.value):
            ptr = self.fabm.get_variable(self.pmodel, SCALAR_DEPENDENCY, i + 1)
            self.scalar_dependencies._data.append(
                Dependency(self, ptr, (), self.fabm.link_scalar)
            )
        for i in range(nconserved.value):
            self.fabm.get_variable_metadata(
                self.pmodel,
                CONSERVED_QUANTITY,
                i + 1,
                ATTRIBUTE_LENGTH,
                strname,
                strunits,
                strlong_name,
                strpath,
            )
            self.conserved_quantities._data.append(
                Variable(
                    self,
                    strname.value.decode("ascii"),
                    strunits.value.decode("ascii"),
                    strlong_name.value.decode("ascii"),
                    strpath.value.decode("ascii"),
                )
            )
        for i in range(nparameters.value):
            self.fabm.get_parameter_metadata(
                self.pmodel,
                i + 1,
                ATTRIBUTE_LENGTH,
                strname,
                strunits,
                strlong_name,
                ctypes.byref(typecode),
                ctypes.byref(has_default),
            )
            self.parameters._data.append(
                Parameter(
                    self,
                    strname.value.decode("ascii"),
                    i,
                    type=typecode.value,
                    units=strunits.value.decode("ascii"),
                    long_name=strlong_name.value.decode("ascii"),
                    has_default=has_default.value != 0,
                )
            )

        self.couplings = NamedObjectList(
            [Coupling(self, i) for i in range(ncouplings.value)]
        )

        # Arrays that combine variables from pelagic and boundary domains.
        self.state_variables = (
            self.interior_state_variables
            + self.surface_state_variables
            + self.bottom_state_variables
        )
        self.diagnostic_variables = (
            self.interior_diagnostic_variables + self.horizontal_diagnostic_variables
        )
        self.dependencies = (
            self.interior_dependencies
            + self.horizontal_dependencies
            + self.scalar_dependencies
        )
        self.variables: NamedObjectList[VariableFromPointer] = (
            self.state_variables + self.diagnostic_variables + self.dependencies
        )

        if settings is not None:
            self._restore_state(settings)

        # For backward compatibility
        self.bulk_state_variables = self.interior_state_variables
        self.bulk_diagnostic_variables = self.interior_diagnostic_variables

        self.itime = -1.0

    def getRates(self, t: Optional[float] = None, surface: bool = True, bottom: bool = True):
        """Returns the local rate of change in state variables,
        given the current state and environment.
        """
        assert self.fabm.idepthdim == -1
        if t is None:
            t = self.itime
        sources = np.empty_like(self._state)
        sources_interior = sources[: len(self.interior_state_variables), ...]
        sources_surface = sources[
            len(self.interior_state_variables) : len(self.interior_state_variables)
            + len(self.surface_state_variables),
            ...,
        ]
        sources_bottom = sources[
            len(self.interior_state_variables) + len(self.surface_state_variables) :,
            ...,
        ]
        assert not (
            (surface or bottom) and self._cell_thickness is None
        ), "You must assign model.cell_thickness to use getRates"
        self.fabm.get_sources(
            self.pmodel,
            t,
            sources_interior,
            sources_surface,
            sources_bottom,
            surface,
            bottom,
            self._cell_thickness,
        )
        if hasError():
            raise FABMException(getError())
        return sources

    def get_sources(
        self,
        t: Optional[float] = None,
        out: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if t is None:
            t = self.itime
        if out is None:
            sources_interior = np.empty_like(self._interior_state)
            sources_surface = np.empty_like(self._surface_state)
            sources_bottom = np.empty_like(self._bottom_state)
        else:
            sources_interior, sources_surface, sources_bottom = out
        assert (
            self._cell_thickness is not None
        ), "You must assign model.cell_thickness to use get_sources"
        self.fabm.get_sources(
            self.pmodel,
            t,
            sources_interior,
            sources_surface,
            sources_bottom,
            True,
            True,
            self._cell_thickness,
        )
        if hasError():
            raise FABMException(getError())
        return sources_interior, sources_surface, sources_bottom

    def get_vertical_movement(self, out: Optional[np.ndarray] = None) -> np.ndarray:
        if out is None:
            out = np.empty_like(self._interior_state)
        self.fabm.get_vertical_movement(self.pmodel, out)
        if hasError():
            raise FABMException(getError())
        return out

    def get_conserved_quantities(self, out: Optional[np.ndarray] = None) -> np.ndarray:
        assert (
            self._cell_thickness is not None
        ), "You must assign model.cell_thickness to use get_conserved_quantities"
        if out is None:
            out = np.empty(
                (len(self.conserved_quantities),) + self.horizontal_domain_shape,
                dtype=self.fabm.numpy_dtype,
            )
        self.fabm.get_conserved_quantities(self.pmodel, out, self._cell_thickness)
        if hasError():
            raise FABMException(getError())
        return out

    def check_state(self, repair: bool = False) -> bool:
        valid = self.fabm.check_state(self.pmodel, repair) != 0
        if hasError():
            raise FABMException(getError())
        return valid

    checkState = check_state

    def getJacobian(self, pert: Union[float, np.ndarray, None] = None):
        # Define perturbation per state variable.
        y_pert = np.empty_like(self.state)
        if pert is None:
            pert = 1e-6
        y_pert[:] = pert

        # Compute dy for original state (used as reference for finite
        # differences later on)
        dy_ori = self.getRates()

        # Create memory for Jacobian
        Jac = np.empty((len(self.state), len(self.state)), dtype=self.state.dtype)

        for i in range(len(self.state)):
            # Save original state variable value, create perturbed one.
            y_ori = self.state[i]
            self.state[i] += y_pert[i]

            # Compute dy for perturbed state, compute Jacobian elements using
            # finite difference.
            dy_pert = self.getRates()
            Jac[:, i] = (dy_pert - dy_ori) / y_pert[i]

            # Restore original state variable value.
            self.state[i] = y_ori

        return Jac

    def findParameter(self, name: str, case_insensitive: bool = False):
        return self.parameters.find(name, case_insensitive)

    def findDependency(self, name: str, case_insensitive: bool = False):
        return self.dependencies.find(name, case_insensitive)

    def findStateVariable(self, name: str, case_insensitive: bool = False):
        return self.state_variables.find(name, case_insensitive)

    def findDiagnosticVariable(self, name: str, case_insensitive: bool = False):
        return self.diagnostic_variables.find(name, case_insensitive)

    def findCoupling(self, name: str, case_insensitive: bool = False):
        return self.couplings.find(name, case_insensitive)

    def find_standard_variable(self, name: str) -> Optional[StandardVariable]:
        pointer = self.fabm.find_standard_variable(name.encode("ascii"))
        if pointer:
            return StandardVariable(self, pointer)

    def require_data(self, standard_variable: StandardVariable):
        return self.fabm.require_data(self.pmodel, standard_variable._pvariable)

    def _get_parameter_tree(self) -> Mapping:
        root = {}
        for parameter in self.parameters:
            pathcomps = parameter.name.split("/")
            parent = root
            for component in pathcomps[:-1]:
                parent = root.setdefault(component, {})
            parent[pathcomps[-1]] = parameter
        return root

    def start(self, verbose: bool = True, stop: bool = False) -> bool:
        ready = True
        if self.fabm.mask_type and self._mask is None:
            log("Mask not yet assigned")
            ready = False
        if self.fabm.variable_bottom_index and self._bottom_index is None:
            log("Bottom indices not yet assigned")
            ready = False

        def process_dependencies(dependencies: Sequence[Dependency]):
            ready = True
            for dependency in dependencies:
                if dependency.required and dependency.value is None:
                    log(f"Value for dependency {dependency.name} is not set.")
                    ready = False
            return ready

        ready = process_dependencies(self.interior_dependencies) and ready
        ready = process_dependencies(self.horizontal_dependencies) and ready
        ready = process_dependencies(self.scalar_dependencies) and ready
        assert ready or not stop, "Not all dependencies have been fulfilled."

        self.fabm.start(self.pmodel)
        if hasError():
            return False
        for i, variable in enumerate(self.interior_diagnostic_variables):
            pdata = self.fabm.get_interior_diagnostic_data(self.pmodel, i + 1)
            if pdata:
                arr = np.ctypeslib.as_array(pdata, self.interior_domain_shape)
                variable._data = arr.view(dtype=self.fabm.numpy_dtype)
            else:
                variable._data = None
        for i, variable in enumerate(self.horizontal_diagnostic_variables):
            pdata = self.fabm.get_horizontal_diagnostic_data(self.pmodel, i + 1)
            if pdata:
                arr = np.ctypeslib.as_array(pdata, self.horizontal_domain_shape)
                variable._data = arr.view(dtype=self.fabm.numpy_dtype)
            else:
                variable._data = None
        return ready

    checkReady = start

    def updateTime(self, nsec: float):
        self.itime = nsec

    def printInformation(self):
        """Show information about the model."""

        def printArray(classname: str, array: Sequence[Variable]):
            if not array:
                return
            log(f" {len(array)} {classname}:")
            for variable in array:
                log(f"    {variable.name} = {variable.value} {variable.units}")

        def parameter2str(p: Parameter) -> str:
            return f"{p.value} {p.units}"

        log("FABM model contains the following:")
        printArray("interior state variables", self.interior_state_variables)
        printArray("bottom state variables", self.bottom_state_variables)
        printArray("surface state variables", self.surface_state_variables)
        printArray("interior diagnostic variables", self.interior_diagnostic_variables)
        printArray(
            "horizontal diagnostic variables", self.horizontal_diagnostic_variables
        )
        printArray("external variables", self.dependencies)
        log(f" {len(self.parameters)} parameters:")
        printTree(self._get_parameter_tree(), parameter2str, "    ")


class Simulator(object):
    def __init__(self, model: Model):
        assert (
            model._cell_thickness is not None
        ), "You must assign model.cell_thickness to use Simulator"
        self.model = model

    def integrate(
        self,
        y0: np.ndarray,
        t: np.ndarray,
        dt: float,
        surface: bool = True,
        bottom: bool = True,
    ):
        y = np.empty((t.size, self.model.state.size))
        self.model.fabm.integrate(
            self.model.pmodel,
            t.size,
            self.model.state.size,
            t,
            y0,
            y,
            dt,
            surface,
            bottom,
            ctypes.byref(self.model._cell_thickness),
        )
        if hasError():
            raise FABMException(getError())
        return y


def unload():
    global ctypes

    for lib in name2lib.values():
        handle = lib._handle
        if os.name == "nt":
            import ctypes.wintypes

            ctypes.windll.kernel32.FreeLibrary.argtypes = [ctypes.wintypes.HMODULE]
            ctypes.windll.kernel32.FreeLibrary(handle)
        else:
            dlclose = ctypes.CDLL(None).dlclose
            dlclose.argtypes = [ctypes.c_void_p]
            dlclose.restype = ctypes.c_int
            dlclose(handle)
    name2lib.clear()


def get_version() -> Optional[str]:
    for lib in name2lib.values():
        version_length = 256
        strversion = ctypes.create_string_buffer(version_length)
        lib.get_version(version_length, strversion)
        return strversion.value.decode("ascii")
