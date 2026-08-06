"""Streamlit components SABES needs and Streamlit does not provide."""
# NOTE: the drawing function is `svg_canvas`, not `canvas`. Binding it to
# `canvas` here would shadow the submodule of the same name, so
# `from sabes.components import canvas` would hand back a function.
from .canvas import (KIND_OPTIC, KIND_POINT, beam, circle, make_spec, rect,
                     svg_canvas)

__all__ = ["svg_canvas", "make_spec", "rect", "circle", "beam",
           "KIND_OPTIC", "KIND_POINT"]
