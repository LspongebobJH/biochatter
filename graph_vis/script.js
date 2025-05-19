const { Graph } = G6;
const graph_input = document.getElementById("graph_input");
const graph_reader = new FileReader();
const render_button = document.getElementById("render_button");
const addNodeButton = document.getElementById("addNodeButton");
const addAPIInput = document.getElementById("addAPIInput");
const productsInput = document.getElementById("productsInput");
const removeNodeButton = document.getElementById("removeNodeButton");
const removeAPIInput = document.getElementById("removeAPIInput");
const addEdgeButton = document.getElementById("addEdgeButton");
const removeEdgeButton = document.getElementById("removeEdgeButton");
const dependenciesInput = document.getElementById("dependenciesInput");
const addSourceInput = document.getElementById("addSourceInput");
const addTargetInput = document.getElementById("addTargetInput");
const argsInput = document.getElementById("argsInput");
const argTypesInput = document.getElementById("argTypesInput");
const removeSourceInput = document.getElementById("removeSourceInput");
const removeTargetInput = document.getElementById("removeTargetInput");
const renderGraphButton = document.getElementById("renderGraphButton");
const downloadGraphButton = document.getElementById("downloadGraphButton");

const dataNames = ['data'];
const nodeKeys = ['api', 'products', 'id', '_deprecated', '_comment'];
const edgeKeys = ['dependencies', 'source', 'target', 'args', 'arg_types', '_deprecated', '_comment'];

let graph = null; // Add graph variable to global scope

render_button.disabled = true;
addNodeButton.disabled = true;
removeNodeButton.disabled = true;
addEdgeButton.disabled = true;
removeEdgeButton.disabled = true;

let graphFileUploaded = false;
let graphData = null;

const enableButtonIfReady = () => {
    if (graphFileUploaded) {
        render_button.disabled = false;
    }
};

// Function to check if both inputs are filled
const updateNodeButtonState = () => {
    addNodeButton.disabled = !(addAPIInput.value.trim());
};

// Function to update remove button state
const updateRemoveButtonState = () => {
    removeNodeButton.disabled = !removeAPIInput.value.trim();
};

// Function to check if all edge inputs are filled
const updateAddEdgeButtonState = () => {
    const allFilled = dependenciesInput.value.trim() &&
                     addSourceInput.value.trim() &&
                     addTargetInput.value.trim() &&
                     argsInput.value.trim() &&
                     argTypesInput.value.trim();
    addEdgeButton.disabled = !allFilled;
};

// Function to check if remove edge inputs are filled
const updateRemoveEdgeButtonState = () => {
    const allFilled = removeSourceInput.value.trim() &&
                     removeTargetInput.value.trim();
    removeEdgeButton.disabled = !allFilled;
};

// Add input event listeners to both fields
addAPIInput.addEventListener('input', updateNodeButtonState);
productsInput.addEventListener('input', updateNodeButtonState);
removeAPIInput.addEventListener('input', updateRemoveButtonState);
dependenciesInput.addEventListener('input', updateAddEdgeButtonState);
addSourceInput.addEventListener('input', updateAddEdgeButtonState);
addTargetInput.addEventListener('input', updateAddEdgeButtonState);
argsInput.addEventListener('input', updateAddEdgeButtonState);
argTypesInput.addEventListener('input', updateAddEdgeButtonState);
removeSourceInput.addEventListener('input', updateRemoveEdgeButtonState);
removeTargetInput.addEventListener('input', updateRemoveEdgeButtonState);

// Add click handler for addNodeButton
addNodeButton.addEventListener('click', () => {
    if (!graph) return; // Ensure graph exists
    
    const apiValue = addAPIInput.value.trim();
    const productsValue = productsInput.value.trim();
    const productsList = productsValue ? productsValue.split('\n').map(p => p.trim()).filter(p => p) : [];
    
    // Add new node to the graph
    graph.addNodeData([
        {
            id: apiValue,
            api: apiValue,
            products: productsList,
            style: {
                size: 10,
                icon: true,
                iconText: apiValue,
                iconFontSize: 5,
                iconFill: '#000',
            }
        }
    ]);
});

// Add click handler for removeNodeButton
removeNodeButton.addEventListener('click', () => {
    if (!graph) return; // Ensure graph exists
    
    const nodeId = removeAPIInput.value.trim();
    
    // Remove the node from the graph
    graph.removeNodeData([nodeId]);
});

// Add click handler for addEdgeButton
addEdgeButton.addEventListener('click', () => {
    if (!graph) return;
    
    const args = JSON.parse(argsInput.value.trim());
    const argTypes = JSON.parse(argTypesInput.value.trim());
    let lineDash = 0;
    const dependenciesValue = dependenciesInput.value.trim();
    const dependenciesList = dependenciesValue ? dependenciesValue.split('\n').map(d => d.trim()).filter(d => d) : [];

    const argName = Object.values(args)[0]; // fow now only one edge required arg being considered 
    if (!dataNames.includes(argName)) {
        lineDash = 1;
    };

    graph.addEdgeData([{
        source: addSourceInput.value.trim(),
        target: addTargetInput.value.trim(),
        dependencies: dependenciesList,
        args: args,
        arg_types: argTypes,
        style: {
            lineDash: lineDash,
        }
    }]);
});

// Add click handler for removeEdgeButton
removeEdgeButton.addEventListener('click', () => {
    if (!graph) return;
    
    const source = removeSourceInput.value.trim();
    const target = removeTargetInput.value.trim();
    
    // Find and remove the edge
    const edges = graph.getData().edges;
    const edgeToRemove = edges.find(edge => 
        edge.source === source && edge.target === target
    );
    
    if (edgeToRemove) {
        graph.removeEdgeData([edgeToRemove.id]);
    }
});

// Add click handler for renderGraphButton
renderGraphButton.addEventListener('click', () => {
    if (graph) {
        graph.render();
    }
});

// Add click handler for downloadGraphButton
downloadGraphButton.addEventListener('click', () => {
    if (!graph) return;

    // Get current graph data
    const currentNodes = graph.getNodeData();
    const currentEdges = graph.getEdgeData();

    // Remove style key from nodes and edges
    const cleanNodes = currentNodes.map(node => {
        const filteredNode = {};
        nodeKeys.forEach(key => {
            if (key in node) {
                filteredNode[key] = node[key];
            }
        });
        return filteredNode;
    });
    const cleanEdges = currentEdges.map(edge => {
        const filteredEdge = {};
        edgeKeys.forEach(key => {
            if (key in edge) {
                filteredEdge[key] = edge[key];
            }
        });
        return filteredEdge;
    });

    // Create new graph data object
    const newGraphData = {
        ...graphData,
        nodes: cleanNodes,
        edges: cleanEdges
    };

    // Convert to JSON string
    const jsonString = JSON.stringify(newGraphData, null, 2);

    // Create blob and download link
    const blob = new Blob([jsonString], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = 'graph_data.json';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
});

graph_input.onchange = (event) => {
    const file = event.target.files[0];
    if (file) {
        graph_reader.readAsText(file);
    }
};

graph_reader.onload = (event) => {
    graphFileUploaded = true;
    graphData = JSON.parse(event.target.result);
    enableButtonIfReady();
};

render_button.onclick = () => {
    if (graphData) {
        const nodes = [];
        const edges = [];

        for (let i = 0; i < graphData.nodes.length; i++) {
            if (!("_deprecated" in Object.keys(graphData.nodes[i])) || graphData.nodes[i]._deprecated == false) {
                graphData.nodes[i].style = {
                    size: 10,
                    icon: true,
                    iconText: graphData.nodes[i].id,
                    iconFontSize: 5,
                    iconFill: '#000',
                };
                const node = graphData.nodes[i];
                nodes.push(node);
            }
        }
        
        for (let i = 0; i < graphData.edges.length; i++) {
            if (!("_deprecated" in Object.keys(graphData.edges[i])) || graphData.edges[i]._deprecated == false) {
                const edgeArg = graphData.edges[i].args;
                const argName = Object.values(edgeArg)[0]; // fow now only one edge required arg being considered 
                if (!dataNames.includes(argName)) {
                    graphData.edges[i].style = {
                        lineDash: 1,
                    };
                };
                const edge = graphData.edges[i];
                edges.push(edge);
            }
        }

        let processedGraphData = {};
        Object.assign(processedGraphData, graphData);
        processedGraphData.nodes = nodes;
        processedGraphData.edges = edges;
        console.log(processedGraphData);
        showGraph(processedGraphData);
    }
};

var showGraph = data => {
    graph = new Graph({
        container: 'container',
        autoFit: 'view',
        autoResize: true,
        data,
        node: {
            palette: {
                field: 'group',
                color: 'tableau',
            },
        },
        edge: {
            style: {
                haloLineWidth: 8,
                endArrow: true,
                endArrowSize: 2,
            },
        },
        layout: {
            type: 'dendrogram',
            direction: 'LR',
            nodeSep: 20,
            rankSep: 50,
            radial: false,
        },
        behaviors: [
            'zoom-canvas', 'drag-element',
            'hover-activate',
            {
                type: 'click-select',
                key: 'click-select-1',
                multiple: true,
            },
            {
                type: 'drag-canvas',
                key: 'drag-canvas-1',   
            },
        ],
        plugins: [
            {
                type: 'tooltip',
                trigger: 'hover',
                getContent: (e, _item) => {
                    const { targetType, target } = e;
                    const item = _item[0];
                    if (targetType == 'node') {
                        let content = `
                            <div style="word-wrap:break-word;">
                                <p>API: ${item.api}<p>
                                <p>products: ${item.products}<p>
                            </div>
                        `;
                        return content;
                    } else if (targetType == 'edge') {
                        let content = `
                            <div style="word-wrap:break-word;">
                                <p>dependencies: ${item.dependencies}<p>
                                <p>args: ${JSON.stringify(item.args)}<p>
                                <p>arg_types: ${JSON.stringify(item.arg_types)}<p>
                            </div>
                        `;
                        return content;
                    }
                    return null;
                }
            },
            {
                key: 'minimap',
                type: 'minimap',
                size: [240, 160],
                position: 'right-bottom'
            },
        ]
    });
    
    graph.render();
}; 