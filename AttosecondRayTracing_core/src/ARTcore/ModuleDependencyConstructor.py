import networkx as nx
import logging
import numpy as np


# Configure logging
# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

class DepSolver:
    def __init__(self, nodes, functions):
        self.G = nx.DiGraph()
        self.nodes = nodes
        self.functions = functions

        # Add nodes to the graph
        for node in nodes:
            self.G.add_node(node)

        # Add edges and functions to the graph
        for outputs, func_data in functions.items():
            for inputs, func in func_data:
                for output in outputs:
                    for input_node in inputs:
                        if self.G.has_edge(input_node, output):
                            self.G.edges[input_node, output]['functions'].append((func, inputs, outputs))
                        else:
                            self.G.add_edge(input_node, output, functions=[(func, inputs, outputs)])

    def calculate_values(self, verify_consistency=False, **known_values):
        known_values = {k: v for k, v in known_values.items() if v is not None}
        G_copy = self.G.copy()
        steps_to_find_value = {node: None for node in G_copy.nodes}

        # Set known values in the graph
        for node, value in known_values.items():
            if node in G_copy.nodes:
                G_copy.nodes[node]['value'] = value
                steps_to_find_value[node] = 0

        iterations = 0
        max_iterations = len(G_copy.nodes) * 10  # Prevent infinite loops

        while iterations < max_iterations:
            calculated = False
            for node in G_copy.nodes:
                if 'value' not in G_copy.nodes[node]:
                    iterations += 1
                    for pred in G_copy.predecessors(node):
                        edge_data = G_copy.edges[pred, node]
                        for func, inputs, outputs in edge_data['functions']:
                            if all('value' in G_copy.nodes[inp] for inp in inputs):
                                inputs_values = [G_copy.nodes[inp]['value'] for inp in inputs]
                                try:
                                    results = func(*inputs_values)
                                except Exception as e:
                                    logger.error(f"Error while calculating {node} using {func}: {e}")
                                    raise ValueError(f"Error while calculating {node} using {func}: {e}")

                                if not isinstance(results, tuple):
                                    results = (results,)
                                for output, result in zip(outputs, results):
                                    if 'value' in G_copy.nodes[output] and not self._is_close(G_copy.nodes[output]['value'], result):
                                        conflict_message = f"Conflict detected: {output} cannot be both {G_copy.nodes[output]['value']} and {result}."
                                        logger.error(conflict_message)
                                        raise ValueError(conflict_message)
                                    G_copy.nodes[output]['value'] = result
                                    steps_to_find_value[output] = iterations + 1
                                    calculated = True
                                    logger.debug(f"Calculated {output} as {result} using {func.__name__} with inputs {inputs_values}")

            if not calculated:
                break

        if iterations == max_iterations:
            unresolved_nodes = [node for node in G_copy.nodes if 'value' not in G_copy.nodes[node]]
            error_message = f"Max iterations reached. The graph might have a cycle or unresolved dependencies. Unresolved nodes: {unresolved_nodes}"
            logger.error(error_message)
            raise RuntimeError(error_message)

        # Check for unresolved nodes
        unresolved_nodes = [node for node in G_copy.nodes if 'value' not in G_copy.nodes[node] or G_copy.nodes[node]['value'] is None]
        if unresolved_nodes:
            error_message = f"Not enough data to fill in the graph. Unresolved nodes: {unresolved_nodes}"
            logger.error(error_message)
            raise ValueError(error_message)

        # Log final state of the graph
        logger.debug(f"Final graph state: {G_copy.nodes(data=True)}")

        # Verify consistency of the calculated values
        if verify_consistency:
            self.verify_consistency(G_copy)

        # Return the calculated values along with the steps to find each value
        return {node: data['value'] for node, data in G_copy.nodes(data=True) if 'value' in data and data['value'] is not None}, steps_to_find_value

    def _is_close(self, a, b, tol=1e-6):
        try:
            return np.isclose(a, b, atol=tol)
        except TypeError:
            return False

    def verify_consistency(self, G):
        tolerance = 1e-6
        for node in G.nodes:
            if 'value' in G.nodes[node] and G.nodes[node]['value'] is not None:
                for succ in G.successors(node):
                    edge_data = G.edges[node, succ]
                    for func, inputs, outputs in edge_data['functions']:
                        inputs_values = [G.nodes[inp]['value'] for inp in inputs]
                        expected_results = outputs
                        try:
                            results = func(*inputs_values)
                        except Exception as e:
                            logger.error(f"Error while verifying {succ} using {func}: {e}")
                            raise ValueError(f"Error while verifying {succ} using {func}: {e}")

                        if not isinstance(results, tuple):
                            results = (results,)

                        for output, result in zip(expected_results, results):
                            if not self._is_close(G.nodes[output]['value'], result, tol=tolerance):
                                error_message = f"Inconsistent result: {output} is {G.nodes[output]['value']} but expected {result} from {func.__name__}."
                                logger.error(error_message)
                                raise ValueError(error_message)

        logger.debug("Consistency check passed for all calculated values.")