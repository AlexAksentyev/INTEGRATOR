# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_dop853_Drift_22_vf', [dirname(__file__)])
        except ImportError:
            import _dop853_Drift_22_vf
            return _dop853_Drift_22_vf
        if fp is not None:
            try:
                _mod = imp.load_module('_dop853_Drift_22_vf', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _dop853_Drift_22_vf = swig_import_helper()
    del swig_import_helper
else:
    import _dop853_Drift_22_vf
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0



def new_doubleArray(*args):
  return _dop853_Drift_22_vf.new_doubleArray(*args)
new_doubleArray = _dop853_Drift_22_vf.new_doubleArray

def delete_doubleArray(*args):
  return _dop853_Drift_22_vf.delete_doubleArray(*args)
delete_doubleArray = _dop853_Drift_22_vf.delete_doubleArray

def doubleArray_getitem(*args):
  return _dop853_Drift_22_vf.doubleArray_getitem(*args)
doubleArray_getitem = _dop853_Drift_22_vf.doubleArray_getitem

def doubleArray_setitem(*args):
  return _dop853_Drift_22_vf.doubleArray_setitem(*args)
doubleArray_setitem = _dop853_Drift_22_vf.doubleArray_setitem

def new_intArray(*args):
  return _dop853_Drift_22_vf.new_intArray(*args)
new_intArray = _dop853_Drift_22_vf.new_intArray

def delete_intArray(*args):
  return _dop853_Drift_22_vf.delete_intArray(*args)
delete_intArray = _dop853_Drift_22_vf.delete_intArray

def intArray_getitem(*args):
  return _dop853_Drift_22_vf.intArray_getitem(*args)
intArray_getitem = _dop853_Drift_22_vf.intArray_getitem

def intArray_setitem(*args):
  return _dop853_Drift_22_vf.intArray_setitem(*args)
intArray_setitem = _dop853_Drift_22_vf.intArray_setitem

def Integrate(*args):
  return _dop853_Drift_22_vf.Integrate(*args)
Integrate = _dop853_Drift_22_vf.Integrate

def InitBasic(*args):
  return _dop853_Drift_22_vf.InitBasic(*args)
InitBasic = _dop853_Drift_22_vf.InitBasic

def CleanUp():
  return _dop853_Drift_22_vf.CleanUp()
CleanUp = _dop853_Drift_22_vf.CleanUp

def InitInteg(*args):
  return _dop853_Drift_22_vf.InitInteg(*args)
InitInteg = _dop853_Drift_22_vf.InitInteg

def ClearInteg():
  return _dop853_Drift_22_vf.ClearInteg()
ClearInteg = _dop853_Drift_22_vf.ClearInteg

def InitEvents(*args):
  return _dop853_Drift_22_vf.InitEvents(*args)
InitEvents = _dop853_Drift_22_vf.InitEvents

def ClearEvents():
  return _dop853_Drift_22_vf.ClearEvents()
ClearEvents = _dop853_Drift_22_vf.ClearEvents

def InitExtInputs(*args):
  return _dop853_Drift_22_vf.InitExtInputs(*args)
InitExtInputs = _dop853_Drift_22_vf.InitExtInputs

def ClearExtInputs():
  return _dop853_Drift_22_vf.ClearExtInputs()
ClearExtInputs = _dop853_Drift_22_vf.ClearExtInputs

def SetRunParameters(*args):
  return _dop853_Drift_22_vf.SetRunParameters(*args)
SetRunParameters = _dop853_Drift_22_vf.SetRunParameters

def ClearParams():
  return _dop853_Drift_22_vf.ClearParams()
ClearParams = _dop853_Drift_22_vf.ClearParams

def Reset():
  return _dop853_Drift_22_vf.Reset()
Reset = _dop853_Drift_22_vf.Reset

def SetContParameters(*args):
  return _dop853_Drift_22_vf.SetContParameters(*args)
SetContParameters = _dop853_Drift_22_vf.SetContParameters

def Vfield(*args):
  return _dop853_Drift_22_vf.Vfield(*args)
Vfield = _dop853_Drift_22_vf.Vfield

def Jacobian(*args):
  return _dop853_Drift_22_vf.Jacobian(*args)
Jacobian = _dop853_Drift_22_vf.Jacobian

def JacobianP(*args):
  return _dop853_Drift_22_vf.JacobianP(*args)
JacobianP = _dop853_Drift_22_vf.JacobianP

def AuxFunc(*args):
  return _dop853_Drift_22_vf.AuxFunc(*args)
AuxFunc = _dop853_Drift_22_vf.AuxFunc

def MassMatrix(*args):
  return _dop853_Drift_22_vf.MassMatrix(*args)
MassMatrix = _dop853_Drift_22_vf.MassMatrix
# This file is compatible with both classic and new-style classes.


