#!/usr/bin/env python3
import argparse
import json


def main(args):
    cell_w, cell_d = None, None
    if args.cell_dimensions is not None:
        cell_w, cell_d = args.cell_dimensions.split(",")
        cell_w, cell_d = float(cell_w), float(cell_d)

    with open(args.gltf_path) as f:
        gltf = json.load(f)

    src_node_idx, src_node = None, None
    for i, node in enumerate(gltf["nodes"]):
        if node["name"] == args.source_node:
            src_node_idx = i
            src_node = node
            break
    assert src_node is not None
    assert len(src_node.get("children", [])) == 0

    src_mesh_idx = src_node["mesh"]
    src_mesh = gltf["meshes"][src_mesh_idx]
    acc = gltf["accessors"][src_mesh["primitives"][0]["attributes"]["POSITION"]]
    mid_x = (acc["max"][0] - acc["min"][0]) / 2.0
    mid_y = (acc["max"][1] - acc["min"][1]) / 2.0

    if (cell_w, cell_d) == (None, None):
        cell_w = acc["max"][0] - acc["min"][0]
        cell_d = acc["max"][1] - acc["min"][1]
    cell_mid_x = cell_w / 2.0
    cell_mid_y = cell_d / 2.0

    scene = gltf["scenes"][gltf["scene"]]

    scene["nodes"].extend([len(gltf["nodes"]) + i for i in range(len(args.instance))])
    for i, instance in enumerate(args.instance):
        parts = instance.split(",")
        x, y = float(parts[0]), float(parts[1])
        xscale = x * cell_w + cell_mid_x
        yscale = y * cell_d + cell_mid_y
        trans = [xscale / mid_x, yscale / mid_y, 0.0]
        name = f"{args.source_node}_clone{i}" if len(parts) < 3 else parts[2]
        gltf["nodes"].append(
            {
                "name": name,
                "mesh": src_mesh_idx,
                "translation": trans,
                "rotation": src_node.get("rotation", None),
                "scale": src_node.get("scale", None),
            }
        )
    print(json.dumps(gltf, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="clone things")
    parser.add_argument(
        "gltf_path", type=str, help="path to the GLTF scene to use as a base"
    )
    parser.add_argument("source_node", type=str, help="name of the node to clone")
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
