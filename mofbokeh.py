from bokeh.io import output_notebook
from bokeh.models import HoverTool
from bokeh.plotting import figure, show, ColumnDataSource
import pandas as pd
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG
from rdkit.Chem import AllChem

from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*')

# Read the csv file
from ast import literal_eval
df = pd.read_csv("_tables/DataFrameComplete.csv")



df['info.mofid.smiles_linkersList'] = [literal_eval(x) if pd.notna(x) else [] for x in df['info.mofid.smiles_linkers']]
df['info.mofid.smiles_nodesList'] = [literal_eval(x) if pd.notna(x) else [] for x in df['info.mofid.smiles_nodes']]

bg = 'outputs.pbe.bandgap'
lcd = 'info.lcd'
cond = 'sigma1_traceAvlog10'
smilesstrings = 'info.mofid.smiles'
col = 'color'
df = df.dropna(subset = [smilesstrings, bg, lcd, cond])

df = df.sample(frac = 0.3)

#_______________________Colormapper

df['color'] = '#A9A9A9'

colormap = {'Fe': '#EF4444', 'Nd': '#D54799','Co':'#FAA31B','Mn':'#394BA0','Zn': '#009F75', 'Others':'#A9A9A9'}
l = list(colormap.keys())
l.remove('Others')

for elem in l:
    df.loc[df['info.mofid.smiles_nodes'].str.contains(elem), 'color'] = colormap[elem]

colormap_orig = colormap
colormap = {val: key for key, val in colormap.items()} #Reverse the color map

#___________________RDKit molecules

drop = []
imgs = []
for (i, row) in df.iterrows():
    
    smiles = row[smilesstrings]
    try:
        mol = MolFromSmiles(smiles)
        d2d = MolDraw2DSVG(150, 150)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        svg = d2d.GetDrawingText()
        imgs.append(svg)
    except:
        drop.append(i)

df.drop(drop, inplace=True)
df['img'] = imgs

# Get data to plot
# df = df[df['color'].isin(['red', 'purple'])]

all_smiles = df[smilesstrings]
x = df[bg].values
# x = df[lcd].values
y = df[cond].values
color = df[col].values
imgs = df['img']


#__________________Draw outputs

output_notebook(hide_banner=True)

p = figure()

points = {}

for i, dfg in df.groupby('color'):
    all_smiles = dfg[smilesstrings]
    x = dfg[bg].values
    # x = dfg[lcd].values
    y = dfg[cond].values
    color = dfg[col].values
    imgs = dfg['img']

    source = ColumnDataSource(data={"x": x,
                                    "y": y,
                                    "imgs": imgs,
                                    "smiles": all_smiles,
                                    "colors": color
                                    }
                                    )
    metal = colormap[i]
    points[metal] = p.scatter("x", "y", color = 'colors', size =10, source=source, legend_label = colormap[i], alpha = 0.1, muted_alpha=0.9)

    # Create tooltips referencing stored images
    TOOLTIPS = """\
        <div>
            <div>
                @imgs{safe}
            </div>
            <div>
                <span>[$index]</span>
            </div>
            <div>
                <span>@smiles{safe}</span>
            </div>
            <div>
                <span>($x, $y)</span>
            </div>
        </div>
    """

    # Connect tooltips to plot
    p.add_tools(HoverTool(tooltips=TOOLTIPS))



p.height = 600
p.width = 600
#p.sizing_mode = "scale_width"
p.xaxis.axis_label = "Bandgap (eV)"
p.yaxis.axis_label = "log of Conductivity (S/m)"
p.xaxis.axis_label_text_font_size = '18pt'
p.yaxis.axis_label_text_font_size = '18pt'

p.legend.location = "top_left"
p.legend.click_policy="mute"

show(p)
    