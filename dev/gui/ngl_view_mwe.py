from pyhdx.panel.widgets import NGLView_factory
import param
from pathlib import Path
import panel as pn


class watcher(object):
    def __init__(self, ngl_viewer, parent):
        self.ngl_viewer = ngl_viewer
        self.parent = parent
        self.parent.param.watch(self.changeSpin, "spin")
        self.parent.param.watch(self.changeRepresentation, "representation")
        self.parent.param.watch(self.changerscb, "rcsb_id")
        #self.parent.param.watch(self.changeColor, "no_coverage")
        self.parent.param.watch(self.changeColorList, "color_list")
        self.parent.param.watch(self.changePDBString, "pdb_string")


    def changeSpin(self, event):
        self.ngl_viewer.spin = event.new


    def changeRepresentation(self, event):
        print(event.new)
        self.ngl_viewer.representation = event.new


    def changerscb(self, event):
        self.ngl_viewer.rscb = event.new


    def changeColor(self, event):
        self.ngl_viewer.color = event.new.no_coverage.replace('#', '0x')


    def changeColorList(self, event):
        self.ngl_viewer.color_list = event.new

    def changePDBString(self, event):
        print('pdb changed')
        self.ngl_viewer.pdb_string = event.new


class Test(param.Parameterized):
    pdb_string = param.String()
    rcsb_id = param.String(default="1qyn", doc='RCSB PDB identifier of protein entry to download and visualize.')
    representation = param.Selector(default='cartoon',
                                    objects=['backbone', 'ball+stick', 'cartoon', 'hyperball', 'licorice',
                                             'ribbon', 'rope', 'spacefill', 'surface'],
                                    doc='Representation to use to render the protein.')
    spin = param.Boolean(default=False, doc='Rotate the protein around an axis.')
    #no_coverage = param.Color(default='#8c8c8c', doc='Color to use for regions of no coverage.')
    color_list = param.List(default=[["red", "64-74 or 134-154 or 222-254 or 310-310 or 322-326"],
                                     ["green", "311-322"],
                                     ["yellow",
                                      "40-63 or 75-95 or 112-133 or 155-173 or 202-221 or 255-277 or 289-309"],
                                     ["blue", "1-39 or 96-112 or 174-201 or 278-288"],
                                     ["white", "*"]])


pdb_string = Path('1qyn.pdb').read_text()


view = NGLView_factory.create_view(pdb_string=pdb_string, sizing_mode='stretch_both')
t = Test()
watch = watcher(view, t)

app = pn.Row(t.param, view)
app.servable()


def cb():
    view.pdb_string = pdb_string

pn.state.onload(cb)


#result.servable()