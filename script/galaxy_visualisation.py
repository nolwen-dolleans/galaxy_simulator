import numpy as np
from vispy import app, scene
from vispy.color import get_colormap
import pathlib
import sys
import os
from tqdm import tqdm


def load_data(file_path = "results/Barnes_Hut_parallelstars_data"):
    
    """
    size = os.path.getsize(file_path)
    with open(file_path, "rb") as f, tqdm(total=size, unit="B", unit_scale=True, desc="Load Data: ") as bar:
        header = f.readline()
        n=int(header)
        bar.update(len(header))
        data = []
        for line in f:
            data.append(np.fromstring(line, dtype=np.float32, sep = ",")[:3 * n])
            bar.update(len(line))
    
    results = np.stack(data).reshape(-1, n, 3)
    
    np.save("data.npy", results)
    """
    raw = np.load(file_path, mmap_mode='r')   # (T, 3, N) : blocs x, y, z écrits par le C++
    results = raw.transpose(0, 2, 1)          # (T, N, 3) : vue, sans copie
    print(results.shape)
    return results



path = "results/data.npy"
data = load_data(path)


canvas = scene.SceneCanvas(keys="interactive", size=(900, 700),
                        bgcolor="black", show=True)
view = canvas.central_widget.add_view()
view.camera = "turntable"


stars = scene.visuals.Markers(parent=view.scene)
stars.set_data(data[0], size=3, edge_width=0, face_color=(1, 1, 1, 0.3))

stars.set_gl_state("translucent", depth_test=False)

view.camera.set_range()

### -- COLORS -- ###
import matplotlib
cmap = matplotlib.colormaps["magma"]                             


i = 0
playing = True

@canvas.events.key_press.connect
def on_key(event):
    global playing
    if event.key == "Space":
        playing = not playing


def on_timer(event):
    global i
    r = np.hypot(data[i, :, 0], data[i, :, 1])   # rayon de chaque étoile à t=0                    
    v = np.clip(r / np.percentile(r, 98), 0, 1)           
    colors = cmap(0.75 + 0.25 * (1-v)).astype(np.float32)                     
    colors[:, 3] = 0.5
    if not playing:
        return
    i = (i + 1) % len(data)
    stars.set_data(data[i], size=3, edge_width=0, face_color=colors)


timer = app.Timer(interval=1 / 60, connect=on_timer, start=True)

if __name__ == "__main__":
    
    app.run()