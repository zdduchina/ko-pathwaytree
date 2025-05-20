import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import json
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch, Arc
import matplotlib.colors as mcolors
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from typing import Dict, List, Tuple, Union
import argparse

def create_tree_plot(gene_data_file: str, pathway_hierarchy_file: str, output_image: str) -> None:
    """
    Create a hierarchical visualization of KEGG pathways with gene regulation information.
    Only shows labels for the outermost (pathway) layer.
    
    Parameters:
    -----------
    gene_data_file : str
        Path to the gene data file containing pathway IDs and regulation information
    pathway_hierarchy_file : str
        Path to the file containing KEGG pathway hierarchy information
    output_image : str
        Path where the output image will be saved
    """
    # Read data files with explicit data types
    try:
        gene_df = pd.read_csv(gene_data_file, sep='\t')
        pathway_df = pd.read_csv(pathway_hierarchy_file, sep='\t', dtype={
            'level1': str,
            'level2': str,
            'pathway_name': str,
            'pathway_id': str
        })
    except Exception as e:
        raise ValueError(f"Error reading input files: {str(e)}")
    
    # Clean column names
    gene_df.columns = gene_df.columns.str.strip()
    pathway_df.columns = pathway_df.columns.str.strip()
    
    # Create directed graph
    G = nx.DiGraph()
    root = "KEGG Pathways"
    G.add_node(root)
    
    # Create hierarchy dictionary
    hierarchy = {
        'level1': {},  # level1 -> set of level2
        'level2': {},  # level2 -> set of pathways
        'pathways': {} # pathway_name -> pathway_id
    }
    
    # Collect hierarchy relationships
    for _, row in pathway_df.iterrows():
        try:
            level1 = str(row['level1']).strip()
            level2 = str(row['level2']).strip()
            pathway_name = f"{str(row['pathway_name']).strip()} ({str(row['pathway_id']).strip()})"
            pathway_id = str(row['pathway_id']).strip()
            
            if level1 not in hierarchy['level1']:
                hierarchy['level1'][level1] = set()
            hierarchy['level1'][level1].add(level2)
            
            if level2 not in hierarchy['level2']:
                hierarchy['level2'][level2] = set()
            hierarchy['level2'][level2].add(pathway_name)
            
            hierarchy['pathways'][pathway_name] = pathway_id
            
        except Exception as e:
            print(f"Warning: Skipping invalid pathway entry: {str(e)}")
            continue
    
    # Build graph structure
    # 1. Add first layer
    for level1 in sorted(hierarchy['level1'].keys()):
        G.add_edge(root, level1)
    
    # 2. Add second layer
    for level1, level2_set in hierarchy['level1'].items():
        for level2 in sorted(level2_set):
            G.add_edge(level1, level2)
    
    # 3. Add pathway layer
    for level2, pathway_set in hierarchy['level2'].items():
        for pathway in sorted(pathway_set):
            G.add_edge(level2, pathway)
            G.nodes[pathway]['pathway_id'] = hierarchy['pathways'][pathway]
    
    # Calculate gene counts for each pathway
    pathway_counts: Dict[str, Dict[str, int]] = {}
    for _, row in gene_df.iterrows():
        if pd.notna(row['pathway_ids']) and row['pathway_ids'] != '-':
            pathways = str(row['pathway_ids']).split(',')
            for pathway_id in pathways:
                pathway_id = str(pathway_id).strip()
                if pathway_id not in pathway_counts:
                    pathway_counts[pathway_id] = {'up': 0, 'down': 0, 'none': 0}
                
                regulation = str(row['up_down']).strip().lower()
                if regulation in ['up', 'down', 'none']:
                    pathway_counts[pathway_id][regulation] += 1
    
    # Create figure
    plt.figure(figsize=(20, 20))
    ax = plt.gca()
    
    # Calculate node depths
    depths = nx.shortest_path_length(G, root)
    
    # Group nodes by depth
    nodes_by_depth: Dict[int, List[str]] = {}
    for node, depth in depths.items():
        if depth not in nodes_by_depth:
            nodes_by_depth[depth] = []
        nodes_by_depth[depth].append(node)
    
    # Create custom layout
    pos: Dict[str, np.ndarray] = {}
    pos[root] = np.array([0.0, 0.0])
    
    # Create node hierarchy tree
    hierarchy_tree = {
        'root': {
            'node': root,
            'children': {},
            'angle_range': (0, 2 * np.pi)
        }
    }
    
    # Build hierarchy tree
    def build_hierarchy_tree():
        level1_nodes = sorted(nodes_by_depth.get(1, []))
        angle_step_l1 = 2 * np.pi / len(level1_nodes)
        
        level2_info = []
        total_level2_nodes = 0
        for l1_node in level1_nodes:
            level2_children = sorted([n for n in G.neighbors(l1_node) 
                                   if depths.get(n) == 2])
            total_level2_nodes += len(level2_children)
            level2_info.extend([(l2, l1_node) for l2 in level2_children])
        
        angle_step_l2 = 2 * np.pi / total_level2_nodes
        
        for i, l1_node in enumerate(level1_nodes):
            start_angle = i * angle_step_l1
            end_angle = (i + 1) * angle_step_l1
            
            hierarchy_tree['root']['children'][l1_node] = {
                'node': l1_node,
                'children': {},
                'angle_range': (start_angle, end_angle)
            }
        
        for i, (l2_node, l1_parent) in enumerate(level2_info):
            l2_start = i * angle_step_l2
            l2_end = (i + 1) * angle_step_l2
            
            hierarchy_tree['root']['children'][l1_parent]['children'][l2_node] = {
                'node': l2_node,
                'children': {},
                'angle_range': (l2_start, l2_end)
            }
            
            level3_children = sorted([n for n in G.neighbors(l2_node) 
                                   if depths.get(n) == 3])
            
            if level3_children:
                angle_step_l3 = (l2_end - l2_start) / len(level3_children)
                
                for k, l3_node in enumerate(level3_children):
                    l3_start = l2_start + k * angle_step_l3
                    l3_end = l2_start + (k + 1) * angle_step_l3
                    
                    hierarchy_tree['root']['children'][l1_parent]['children'][l2_node]['children'][l3_node] = {
                        'node': l3_node,
                        'angle_range': (l3_start, l3_end)
                    }
    
    build_hierarchy_tree()
    
    def assign_positions(tree_node, depth, parent_angle=None):
        node = tree_node['node']
        angle_range = tree_node['angle_range']
        
        if depth == 0:
            pos[node] = np.array([0.0, 0.0])
            current_angle = 0
        else:
            current_angle = (angle_range[0] + angle_range[1]) / 2
            
            if depth == 1:
                distance = depth * 0.8
            elif depth == 2:
                distance = depth * 0.9
            else:
                distance = depth * 0.8
            
            x = float(np.cos(current_angle) * distance)
            y = float(np.sin(current_angle) * distance)
            pos[node] = np.array([x, y], dtype=np.float64)
        
        children = sorted(tree_node.get('children', {}).items())
        for _, child in children:
            assign_positions(child, depth + 1, current_angle)
    
    assign_positions(hierarchy_tree['root'], 0)
    
    # Draw edges
    edge_lines = []
    edge_colors = []
    
    for u, v in G.edges():
        try:
            if u in pos and v in pos:
                u_depth = depths[u]
                v_depth = depths[v]
                
                if v_depth - u_depth == 1:
                    start = pos[u]
                    end = pos[v]
                    
                    if v_depth < 3:
                        mid_point = (start + end) / 2
                        angle_diff = abs(np.arctan2(end[1], end[0]) - np.arctan2(start[1], start[0]))
                        
                        if v_depth == 2:
                            curve_strength = 0.1 * min(1, angle_diff / np.pi)
                        else:
                            curve_strength = 0.15 * (1 + v_depth * 0.2) * min(1, angle_diff / np.pi)
                        
                        normal = np.array([-mid_point[1], mid_point[0]])
                        normal = normal / np.linalg.norm(normal) if np.linalg.norm(normal) > 0 else np.array([0, 0])
                        control_point = mid_point + curve_strength * normal
                        
                        t = np.linspace(0, 1, 20)
                        curve_points = np.array([
                            (1-t)**2 * start + 2*(1-t)*t * control_point + t**2 * end
                            for t in t
                        ])
                        
                        for i in range(len(curve_points)-1):
                            edge_lines.append([curve_points[i], curve_points[i+1]])
                            alpha = 0.2 if v_depth == 2 else max(0.15, 0.4 - v_depth * 0.1)
                            edge_colors.append(mcolors.to_rgba('lightgray', alpha=alpha))
                    
                    else:
                        angle = np.arctan2(end[1], end[0])
                        if angle < 0:
                            angle += 2 * np.pi
                        
                        outer_radius = 2.0
                        arc_center = np.array([
                            outer_radius * np.cos(angle),
                            outer_radius * np.sin(angle)
                        ])
                        
                        edge_lines.append([start, arc_center])
                        edge_colors.append(mcolors.to_rgba('lightgray', alpha=0.3))
        
        except Exception as e:
            print(f"Warning: Could not draw edge between {u} and {v}: {str(e)}")
            continue

    if edge_lines:
        edge_collection = LineCollection(edge_lines, colors=edge_colors, 
                                      linewidths=0.5,
                                      linestyles='-',
                                      zorder=1)
        ax.add_collection(edge_collection)
    
    # Set color mapping
    level_colors = {
        1: ['#4A90E2', '#8E44AD', '#2ECC71', '#34495E', '#E67E22'],
        2: ['#7FB3D5', '#BB8FCE', '#82E0AA', '#85929E', '#F5B041']
    }
    
    def create_color_gradient(up_ratio: float, down_ratio: float, 
                            none_ratio: float, num_points: int = 50) -> np.ndarray:
        colors = np.zeros((num_points-1, 4))
        for i in range(num_points-1):
            pos_val = i / (num_points - 1)
            if pos_val < up_ratio:
                colors[i] = mcolors.to_rgba('#E41A1C', alpha=0.8)
            elif pos_val < up_ratio + down_ratio:
                colors[i] = mcolors.to_rgba('#377EB8', alpha=0.8)
            else:
                colors[i] = mcolors.to_rgba('#D3D3D3', alpha=0.8)
        return colors

    # Assign colors to nodes
    node_colors: Dict[str, str] = {}
    for depth in [1, 2]:
        if depth in nodes_by_depth:
            nodes = sorted(nodes_by_depth[depth])
            colors = level_colors[depth]
            for i, node in enumerate(nodes):
                node_colors[node] = colors[i % len(colors)]
    
    # Draw nodes and labels
    for node, (x, y) in pos.items():
        try:
            depth = depths[node]
            
            if depth == 1:  # First layer (inner)
                angle = np.arctan2(y, x)
                if angle < 0:
                    angle += 2 * np.pi
                
                # Add label for first layer
                label_radius = 0.9  # Inner radius for level 1
                label_x = label_radius * np.cos(angle)
                label_y = label_radius * np.sin(angle)
                
                # Determine text alignment based on angle
                if -np.pi/2 <= angle <= np.pi/2:
                    ha = 'left'
                    rotation = np.degrees(angle)
                else:
                    ha = 'right'
                    rotation = np.degrees(angle + np.pi)
                
                plt.text(label_x, label_y, node,
                        rotation=rotation,
                        horizontalalignment=ha,
                        verticalalignment='center',
                        fontsize=8,  # Larger font size for level 1
                        alpha=0.8)  # More visible for level 1
                
                # Draw node circle
                num_neighbors = len(list(nx.neighbors(G, node)))
                base_size = 0.15
                node_size = base_size + num_neighbors * 0.01
                
                circle = Circle((x, y), node_size, 
                              facecolor=node_colors.get(node, '#D3D3D3'),
                              edgecolor='none',
                              alpha=0.9)
                ax.add_patch(circle)
            
            elif depth == 2:  # Second layer (middle)
                angle = np.arctan2(y, x)
                if angle < 0:
                    angle += 2 * np.pi
                
                # Add label for second layer with smaller font
                label_radius = 1.8  # Middle radius for level 2
                label_x = label_radius * np.cos(angle)
                label_y = label_radius * np.sin(angle)
                
                # Determine text alignment based on angle
                if -np.pi/2 <= angle <= np.pi/2:
                    ha = 'left'
                    rotation = np.degrees(angle)
                else:
                    ha = 'right'
                    rotation = np.degrees(angle + np.pi)
                
                plt.text(label_x, label_y, node,
                        rotation=rotation,
                        horizontalalignment=ha,
                        verticalalignment='center',
                        fontsize=6,  # Smaller font size for level 2
                        alpha=0.6)  # More transparent for level 2
                
                # Draw node circle
                num_neighbors = len(list(nx.neighbors(G, node)))
                base_size = 0.08
                node_size = base_size + num_neighbors * 0.01
                
                circle = Circle((x, y), node_size, 
                              facecolor=node_colors.get(node, '#D3D3D3'),
                              edgecolor='none',
                              alpha=0.4)
                ax.add_patch(circle)
            
            elif depth == 3:  # Third layer (outermost) - no labels
                pathway_id = G.nodes[node].get('pathway_id')
                if not pathway_id:
                    continue
                    
                counts = pathway_counts.get(pathway_id, {'up': 0, 'down': 0, 'none': 0})
                total = sum(counts.values())
                
                if total == 0:
                    continue
                
                angle = np.arctan2(y, x)
                if angle < 0:
                    angle += 2 * np.pi
                
                inner_radius = 1.45
                outer_radius = 2.0
                
                start_x, start_y = x, y
                end_x = outer_radius * np.cos(angle)
                end_y = outer_radius * np.sin(angle)
                
                up_ratio = counts['up'] / total
                down_ratio = counts['down'] / total
                none_ratio = counts['none'] / total
                
                t = np.linspace(0, 1, 50)
                x_coords = start_x + (end_x - start_x) * t
                y_coords = start_y + (end_y - start_y) * t
                points = np.column_stack([x_coords, y_coords])
                
                segments = np.zeros((len(points)-1, 2, 2))
                segments[:, 0] = points[:-1]
                segments[:, 1] = points[1:]
                
                colors = create_color_gradient(up_ratio, down_ratio, none_ratio)
                
                lc = LineCollection(segments, colors=colors, linewidth=3, zorder=3)
                ax.add_collection(lc)
                
                arc = Arc((0, 0), 
                         outer_radius * 2, outer_radius * 2,
                         theta1=np.degrees(angle - 0.15),
                         theta2=np.degrees(angle + 0.15),
                         color='#666666',
                         linewidth=6,
                         alpha=0.3)
                ax.add_patch(arc)
                
                # No labels for the third layer
        
        except Exception as e:
            print(f"Warning: Could not draw node {node}: {str(e)}")
            continue
    
    # Set figure range and style
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    ax.set_aspect('equal')
    plt.axis('off')
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='Level 1',
                  markerfacecolor=level_colors[1][0], markersize=12, markeredgewidth=0),
        plt.Line2D([0], [0], marker='o', color='w', label='Level 2',
                  markerfacecolor=level_colors[2][0], markersize=10, markeredgewidth=0),
        plt.Line2D([0], [0], color='#E41A1C', label='Up-regulated', linewidth=3),
        plt.Line2D([0], [0], color='#377EB8', label='Down-regulated', linewidth=3),
        plt.Line2D([0], [0], color='#D3D3D3', alpha=0.3, label='No change', linewidth=1.5)
    ]
    plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.1, 1.1))
    
    try:
        plt.savefig(output_image, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Tree plot saved as: {output_image}")
    except Exception as e:
        raise RuntimeError(f"Error saving plot to {output_image}: {str(e)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a KEGG pathway tree visualization with outer labels only')
    parser.add_argument('gene_data', help='Path to gene data file')
    parser.add_argument('pathway_hierarchy', help='Path to pathway hierarchy file')
    parser.add_argument('output', help='Path for output image')
    
    args = parser.parse_args()
    
    try:
        create_tree_plot(args.gene_data, args.pathway_hierarchy, args.output)
    except Exception as e:
        print(f"Error: {str(e)}")
        exit(1) 