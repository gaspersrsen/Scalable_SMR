import openmc
import numpy as np

def filled(x,y,r):
    return x**2 + y**2 - r**2 <= r-1

def fatfilled(x,y,radius):
	return filled(x, y, radius) and not (
		filled(x + 1, y, radius) and
		filled(x - 1, y, radius) and
		filled(x, y + 1, radius) and
		filled(x, y - 1, radius) and
		filled(x + 1, y + 1, radius) and
		filled(x + 1, y - 1, radius) and
		filled(x - 1, y - 1, radius) and
		filled(x - 1, y + 1, radius))

def get_n_fe(d):
    n_instances = 0
    r = (d-1) / 2
    for y in np.linspace(-r,r,d):
        for x in np.linspace(-r,r,d):
            if x >= 0 and y >= 0 and y <= x:
                if filled(x,y,r):
                    n_fei = 8
                    if x==0: n_fei/=2
                    if y==0: n_fei/=2
                    if x==y: n_fei/=2
                    n_instances += int(n_fei)
    return n_instances

def n_neibhours(x,y,radius):
	return (
		int(filled(x + 1, y, radius)) +
		int(filled(x - 1, y, radius)) +
		int(filled(x, y + 1, radius)) +
		int(filled(x, y - 1, radius)))
 
def generate_octant_subchannels_and_connections(d):
    """
    Dynamically maps subchannels and generates lateral crossflow connections
    for the N-NE octant based on grid dimension d.
    """
    r = (d - 1) / 2
    
    # We need exact integer index steps to map a discrete subchannel layout.
    # Convert linspace limits to structured loop coordinates matching your scale
    x_range = np.arange(-r, r + 1)
    y_range = np.arange(-r, r + 1)
    
    # Dictionary to map continuous (x, y) coordinates to a sequential subchannel ID
    subchannel_map = {}
    subchannel_list = []
    subchannel_id = 0
    
    # ==========================================================
    # STEP 1: DISCOVER AND CATALOG ALL N-NE OCTANT SUBCHANNELS
    # ==========================================================
    # Loop over Y first (rows) then X (columns) to match your bottom-to-top layout
    for y in y_range:
        for x in x_range:
            # Enforce your octant boundary constraints
            if x >= 0 and y >= 0 and y <= x:
                if filled(x, y, r):
                    subchannel_map[(x, y)] = subchannel_id
                    subchannel_list.append((x, y))
                    subchannel_id += 1

    # ==========================================================
    # STEP 2: DYNAMICALLY FIND NEIGHBORS AND BUILD CONNECTIONS
    # ==========================================================
    subchannel_connections = []
    
    for (x, y), current_id in subchannel_map.items():
        # Check right neighbor (Horizontal Connection)
        right_neighbor = (x + 1, y)
        if right_neighbor in subchannel_map:
            neighbor_id = subchannel_map[right_neighbor]
            subchannel_connections.append((current_id, neighbor_id))
            
        # Check upper neighbor (Vertical Connection)
        upper_neighbor = (x, y + 1)
        if upper_neighbor in subchannel_map:
            neighbor_id = subchannel_map[upper_neighbor]
            subchannel_connections.append((current_id, neighbor_id))
            
        # Handle special Octant Symmetry Reflection boundaries:
        # If we are sitting exactly on the diagonal (y == x), the upper neighbor (x, y+1)
        # falls outside this specific N-NE octant slice. However, physically, it is a
        # symmetry reflection copy of its rightward neighbor (x+1, y).
        if x == y and x > 0:
            upper_reflected_neighbor = (x, y + 1)
            # Find the ID of the rightward neighbor within our tracked octant map
            if (x + 1, y) in subchannel_map:
                reflected_id = subchannel_map[(x + 1, y)]
                # Connect the diagonal core boundary to its cross-octant mirror reflection
                subchannel_connections.append((current_id, reflected_id))
	
    return subchannel_list, subchannel_connections
