#!/usr/bin/env python3
import argparse

from pygltflib import GLTF2, Node


def main(args):
    cell_w, cell_d = None, None
    if args.cell_dimensions is not None:
        cell_w, cell_d = args.cell_dimensions.split(",")
        cell_w, cell_d = float(cell_w), float(cell_d)

    gltf = GLTF2().load(args.gltf_path)

    src_mesh_idx, src_mesh = None, None
    for i, mesh in enumerate(gltf.meshes):
        if mesh.name == args.source_mesh:
            src_mesh_idx = i
            src_mesh = mesh
            break
    assert src_mesh is not None
    acc = gltf.accessors[src_mesh.primitives[0].attributes.POSITION]
    mid_x = (acc.max[0] - acc.min[0]) / 2.0
    mid_y = (acc.max[1] - acc.min[1]) / 2.0

    if (cell_w, cell_d) == (None, None):
        cell_w = acc.max[0] - acc.min[0]
        cell_d = acc.max[1] - acc.min[1]
    cell_mid_x = cell_w / 2.0
    cell_mid_y = cell_d / 2.0

    scene = gltf.scenes[gltf.scene]

    scene.nodes.extend([len(gltf.nodes) + i for i in range(len(args.instance))])
    for i, instance in enumerate(args.instance):
        parts = instance.split(",")
        x, y = float(parts[0]), float(parts[1])
        xscale = x * cell_w + cell_mid_x
        yscale = y * cell_d + cell_mid_y
        trans = [xscale / mid_x, yscale / mid_y, 0.0]
        name = f"{args.source_mesh}_clone{i}" if len(parts) < 3 else parts[2]
        gltf.nodes.append(Node(name=name, mesh=src_mesh_idx, translation=trans))
    print(gltf.to_json(indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="clone things")
    parser.add_argument(
        "gltf_path", type=str, help="path to the GLTF scene to use as a base"
    )
    parser.add_argument("source_mesh", type=str, help="name of the mesh to clone")
    parser.add_argument(
        "--instance",
        action="append",
        help="x,y[,name] of a clone of the source node",
    )
    parser.add_argument(
        "--cell-dimensions",
        help="width,depth of the grid to use when placing clones."
        " by default width,depth of the source node is used",
    )

    args = parser.parse_args()
    main(args)
