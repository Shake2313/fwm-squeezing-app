"""Explicit CSV alignment actions fence pending graph gestures."""


def test_explicit_alignment_actions_advance_revision_but_rerenders_do_not():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string('''
from types import SimpleNamespace
from unittest.mock import patch
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from gabes_ui.plotcard import _render_experimental_comparison

x = np.linspace(-5., 5., 101)
y = 1. - .3 * np.exp(-x * x)
csv = "\\n".join(f"{a},{b}" for a, b in zip(x, y)).encode()
uploaded = SimpleNamespace(name='trace.csv', getvalue=lambda: csv)
fig, ax = plt.subplots()
ax.plot(x, y)
ax.set_xlim(-10, 10)
view = dict(figure=fig, comparison=dict(axis_index=0, x_unit='MHz'))
with patch('streamlit.delta_generator.DeltaGenerator.file_uploader', return_value=uploaded):
    st.session_state['comparison'] = _render_experimental_comparison(view, 'sas', st.container())
plt.close(fig)
''').run(timeout=20)
    assert not app.exception
    revision = app.session_state["comparison"]["alignment_revision"]
    app.run(timeout=20)
    assert app.session_state["comparison"]["alignment_revision"] == revision
    for key in ("sas___csv_apply", "sas___csv_reset", "sas___csv_frame"):
        next(button for button in app.button if button.key == key).click().run(timeout=20)
        assert not app.exception
        revision += 1
        assert app.session_state["comparison"]["alignment_revision"] == revision
        app.run(timeout=20)
        assert app.session_state["comparison"]["alignment_revision"] == revision
