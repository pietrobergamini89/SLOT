import holoviews as hv
from bokeh.models import Div
import panel as pn


pn.extension(loading_spinner='dots', loading_color='#00aa41')
hv.extension('bokeh')

title = Div(text="<b>THIS IS A TEST</b>",
                 style={'font-size': '16pt', 'color': 'black'}
            )

pn.Column(title).servable()
