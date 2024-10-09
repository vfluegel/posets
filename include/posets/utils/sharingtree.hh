#pragma once

#include <posets/concepts.hh>

namespace posets::utils {

	template <Vector V>
	class sharingtree {
	private:
		struct st_layer;
		using st_layer_ptr = st_layer*;
		struct st_node;
		using st_node_ptr = st_node*;
		struct st_son;
		using st_son_ptr = st_son*;

		st_layer_ptr firstLayer{ nullptr };
		st_layer_ptr lastLayer{ nullptr };
		st_node_ptr root{ nullptr };
		
		struct st_layer {
			st_node_ptr firstNode{ nullptr };
			st_node_ptr lastNode{ nullptr };

			st_layer_ptr previousLayer{ nullptr };
			st_layer_ptr nextLayer{ nullptr };
		};


		struct st_node {
			st_son_ptr firstSon{ nullptr };
			int val;
			bool isRoot{ false };
			bool isEnd{ false };

			st_node_ptr nextNode{ nullptr };
			// still need the auxiliary values
		};

		struct st_son {
			st_node_ptr node{ nullptr };
			st_son_ptr nextSon{ nullptr };
		};

		/*
		Adds the specified node as child
		*/
		void addSon(st_node_ptr node, st_node_ptr son) {
			st_son_ptr currentSon = node->firstSon;
			// Check if the new node needs to be the first son
			if (son->val > currentSon->node->val) {
				st_son_ptr newSon = new st_son();
				newSon->node = son;
				newSon->nextSon = currentSon;
				node->firstSon = newSon;
			}
			else {
				while (currentSon->nextSon != nullptr && currentSon->node->val > son->val) {
					currentSon = currentSon->nextSon;
				}
				if (currentSon->node->val == son->val) {
					// Constraint on successors is not satisfied -> Throw error or insert children as children of node?
				}
				else {
					st_son_ptr newSon = new st_son();
					newSon->node = son;
					newSon->nextSon = currentSon->nextSon;
					currentSon->nextSon = newSon;
				}
			}
		}

		void removeSon(st_node_ptr node, st_node_ptr son) {
			st_son_ptr currentSon = node->firstSon;
			// Special case: node to remove is the beginning of the successors
			if (son == currentSon->node) {
				node->firstSon = currentSon->nextSon;
			}
			else {
				while (currentSon->nextSon != nullptr && currentSon->nextSon->node->val > son->val) {
					currentSon = currentSon->nextSon;
				}
				if (currentSon->nextSon->node == son) {
					currentSon->nextSon = currentSon->nextSon->nextSon;
				}
			}
			
		}

		st_node_ptr hasSon(st_node_ptr node, int val) {
			st_son_ptr currentSon = node->firstSon;
			while (currentSon->node->val > val) {
				if (currentSon->nextSon == nullptr) return nullptr;
				currentSon = currentSon->nextSon;
			}
			if (currentSon->node->val == val) {
				return currentSon->node;
			}
			else {
				return nullptr;
			}
		}

		bool sameSons(st_node_ptr n1, st_node_ptr n2) {
			st_son_ptr currentSon1 = n1->firstSon;
			st_son_ptr currentSon2 = n2->firstSon;
			bool res = false;
			// Iterate while the node values match
			while (currentSon1->node == currentSon2->node) {
				if (currentSon1->nextSon == nullptr || currentSon2->nextSon == nullptr) {
					// If one of the lists has reached its end, we check whether they are equal (should both be nullptr then)
					// Res only becomes true if both nodes have reached their last son!
					res = currentSon1->nextSon == currentSon2->nextSon;
					break;
				};
				currentSon1 = currentSon1->nextSon;
				currentSon2 = currentSon2->nextSon;
			}
			return res;
		}

		void removeNode(st_layer_ptr layer, st_node_ptr node) {
			st_node_ptr currentNode = layer->firstNode;
			// Special case: node to remove is the beginning of the layer
			if (node == currentNode) {
				layer->firstNode = currentNode->nextNode;
			}
			else {
				while (currentNode->nextNode != nullptr && currentNode->nextNode->val > node->val) {
					currentNode = currentNode->nextNode;
				}
				if (currentNode->nextNode == node) {
					if (currentNode->nextNode == layer->lastNode) {
						layer->lastNode = currentNode;
					}
					else {
						currentNode->nextNode = currentNode->nextNode->nextNode;
					}
				}
			}
		}

		st_node_ptr addNode(st_layer_ptr layer, st_node_ptr node) {
			st_node_ptr currentNode = layer->firstNode;
			//check if new node needs to be the first in the layer
			if (node->val > currentNode->val) {
				node->nextNode = currentNode;
				layer->firstNode = node;
				return node;
			}

			while (currentNode->nextNode != nullptr && currentNode->val > node->val) {
				currentNode = currentNode->nextNode;
			}
			// Check if node with same value and successors as node exists in layer
			while (currentNode->val == node->val) {
				if (sameSons(currentNode, node)) {
					// if it does, the result of the insertion is the existing node
					return currentNode;
				}
				else {
					// else we have to check if there is another node with the same value and do the same check
					currentNode = currentNode->nextNode;
				}
			}
			// if there exists no matching node, we insert
			node->nextNode = currentNode->nextNode;
			currentNode->nextNode = node;
			// Change last node if the node previously was the last
			if (currentNode == layer->lastNode) {
				layer->lastNode = node;
			}
			return node;
		}

		st_layer_ptr addFirstLayer() {
			st_layer_ptr newFirst = new st_layer();
			newFirst->nextLayer = firstLayer;
			firstLayer = newFirst;
			return newFirst;
		}

		st_layer_ptr addLastLayer() {
			st_layer_ptr newLast = new st_layer();
			newLast->previousLayer = lastLayer;
			lastLayer = newLast;
			return newLast;
		}

		void deleteLastLayer() {
			st_layer_ptr newLast = lastLayer->previousLayer;
			newLast->nextLayer = nullptr;
			delete lastLayer;
			lastLayer = newLast;
		}

		// This is a simple standard copy without checks, we might only use the simulation based copy and this method will be obsolete!
		st_node_ptr copy(st_node_ptr node, sharingtree* S, st_layer_ptr destinationLayer) {
			st_node_ptr newNode = new st_node();
			newNode->val = node->val;
			if (node->nextNode != nullptr) {
				st_layer_ptr nextLayer = destinationLayer->nextLayer;
				if (nextLayer == nullptr) {
					nextLayer = S->addLastLayer();
				}
				st_son_ptr son = node->firstSon;
				while (son != nullptr) {
					addSon(newNode, copy(son->node, S, nextLayer));
					son = son->nextSon;
				}
			}
			return addNode(destinationLayer, newNode);
		}

		st_node_ptr copyIfNotSimulated(st_node_ptr node, sharingtree* S, st_layer_ptr destinationLayer) {
			// Check if there is a node that simulates the node we want, if so, just return that node instead
			st_node_ptr checkNode = destinationLayer->firstNode;
			while (checkNode != nullptr && checkNode->val > node->val) {
				if (simulates(checkNode, node)) {
					return checkNode;
				}
				else {
					checkNode = checkNode->nextNode;
				}
			}

			// No simulating node was found, we create a new node and copy its children (if they are not simulated!)
			st_node_ptr newNode = new st_node();
			newNode->val = node->val;
			if (node->nextNode != nullptr) {
				st_layer_ptr nextLayer = destinationLayer->nextLayer;
				if (nextLayer == nullptr) {
					nextLayer = S->addLastLayer();
				}
				st_son_ptr son = node->firstSon;
				while (son != nullptr) {
					addSon(newNode, copyIfNotSimulated(son->node, S, nextLayer));
					son = son->nextSon;
				}
			}
			return addNode(destinationLayer, newNode);
		}


		/*
		* Simulation relation check
		*/

		bool simulates(st_node_ptr n1, st_node_ptr n2) {
			if (n1->firstSon == nullptr) return true; // Not sure about this one, have to figure out EoL representation, but might work becasue only End and Root may have no successors
			// Preliminary check: Condition 1, value must be greater or equal
			if (n1->val < n2->val) return false;
			st_son_ptr s1 = n1->firstSon;
			st_son_ptr s2 = n2->firstSon; //Maybe add a check that they aren't just both empty? Can't be empty bc of ST definition!
			// Check every son of n1
			while (s1 != nullptr) {
				// There was no son of n2 that simulates the son of n1, so we know the node doesn't simulate
				if (s2 == nullptr) return false;

				if (s1->val >= s2->val) {
					// Node is a simulation candidate!
					if (simulates(s1, s2)) {
						// If it simulates, we move on to the next child
						s1 = s1->nextSon;
						// We return to the first child (might not be necessary)
						s2 = n2->firstSon;
					}
					else {
						// Try the next
						s2 = s2->nextSon;
					}
				}
				else {
					// Try the next 
					s2 = s2->nextSon;
				}
			}
			// There was a simulating node found for every son, so the node simulates
			return true;
		}


		/*
		* Union algorithm
		* 
		*/
		sharingtree st_union(sharingtree* S, sharingtree* T) {
			sharingtree U{};
			U.root = new st_node;

			st_son_ptr n_s = S->root->firstSon;
			st_son_ptr n_t = T->root->firstSon;

			if (n_s != nullptr && n_t != nullptr) {
				st_layer_ptr newLayer = U.addLastLayer();
				st_node_ptr newChild{ nullptr };
				while (n_s != nullptr || n_t != nullptr) {
					// Case one: One of the lists is done iterating, copy the nodes without match
					if (n_s == nullptr) {
						newChild = copyIfNotSimulated(n_t->node, &U, newLayer);
						n_t = n_t->nextSon;
					}
					else if (n_t == nullptr) {
						newChild = copyIfNotSimulated(n_s->node, &U, newLayer);
						n_s = n_s->nextSon;
					}
					// Case two: The values are identical, we union the two nodes
					else if (n_s->node->val == n_t->node->val) {
						newChild = node_union(n_s->node, n_t->node, &U, newLayer);
						n_s = n_s->nextSon;
						n_t = n_t->nextSon;
					}
					// Case three: One of the lists is "ahead", we know because the nodes in a layer are ordered
					else if (n_s->node->val > n_t->node->val) {
						newChild = copyIfNotSimulated(n_s->node, &U, newLayer);
						n_s = n_s->nextSon;
					}
					else {
						newChild = copyIfNotSimulated(n_t->node, &U, newLayer);
						n_t = n_t->nextSon;
					}

					//Add as a  son
					if (newChild != nullptr) {
						addSon(U.root, newChild);
					}
				}
			}
			return U;
		}

		st_node_ptr node_union(st_node_ptr n_s, st_node_ptr n_t, sharingtree* U, st_layer_ptr newLayer) {
			st_node_ptr newNode{};
			// Start by inserting EoL if the value is EoL
			if (n_s == nullptr) {
				// TODO
			}
			else {
				newNode = new st_node();
				st_node->val = n_s->val;
				st_layer_ptr nextLayer = newLayer->nextLayer;
				if (nextLayer == nullptr) {
					nextLayer = U->addLastLayer();
				}
				st_son_ptr s_s = n_s->firstSon;
				st_son_ptr s_t = n_t->firstSon;

				while (s_s != nullptr || s_t != nullptr) {
					// Case one: One of the lists is done iterating, copy the nodes without match
					if (s_s == nullptr) {
						newChild = copyIfNotSimulated(s_t->node, U, nextLayer);
						s_t = s_t->nextSon;
					}
					else if (s_t == nullptr) {
						newChild = copyIfNotSimulated(s_s->node, U, nextLayer);
						s_s = s_s->nextSon;
					}
					// Case two: The values are identical, we union the two nodes
					else if (s_s->node->val == s_t->node->val) {
						newChild = node_union(s_s->node, s_t->node, U, nextLayer);
						s_s = s_s->nextSon;
						s_t = s_t->nextSon;
					}
					// Case three: One of the lists is "ahead", we know because the nodes in a layer are ordered
					else if (s_s->node->val > s_t->node->val) {
						newChild = copyIfNotSimulated(s_s->node, U, nextLayer);
						s_s = s_s->nextSon;
					}
					else {
						newChild = copyIfNotSimulated(s_t->node, U, nextLayer);
						s_t = s_t->nextSon;
					}

					//Add as a  son
					if (newChild != nullptr) {
						addSon(newNode, newChild);
					}
				}
				addNode(newLayer, newNode);
			}
			return newNode;
		}


	};
}