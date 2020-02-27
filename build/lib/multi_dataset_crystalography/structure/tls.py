
import numpy

import iotbx.pdb

from scitbx.array_family import flex
from mmtbx.tls.tools import tlso, u_cart_from_tls

def uij_from_tls_vector_and_origin(xyz, tls_vector, origin):
    assert tls_vector.shape == (21,)
    tls_obj = tlso(t=tls_vector[0:6], l=tls_vector[6:12], s=tls_vector[12:21], origin=origin)
    return uij_from_tls_object(xyz=xyz, tls_obj=tls_obj)

def uij_from_tls_object(xyz, tls_obj):
    uij_tls = u_cart_from_tls(sites_cart = flex.vec3_double(xyz),
                              selections = [flex.bool(len(xyz), True)],
                              tlsos      = [tls_obj]                        )
    return numpy.array(uij_tls)

def extract_tls_from_pdb(pdb_file):
    ih = iotbx.pdb.hierarchy.input(pdb_file)
    tls_params = ih.input.extract_tls_params(ih.hierarchy)
    return tls_params

