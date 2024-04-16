// jdh CS224 Spring 2023

import java.awt.*;
import java.util.*;
import java.util.List;

public class Graph {
    List<Node> nodes;

    public Graph() {
        this.nodes = new ArrayList<Node>();
    }

    public void addNode(Node n) {
        this.nodes.add(n);
    }

    public void addEdge(Node n1, Node n2, int capacity) {
        this.addEdge(n1, n2, capacity, 0);
    }

    //----------------------------------------------------------------

    public void addEdge(Node n1, Node n2, int capacity, int flow) {
        Edge e1 = new Edge(n1, n2, capacity, flow);
        assert(flow <= capacity);
        int idx1 = this.nodes.indexOf(n1);
        if (idx1 >= 0) {
            this.nodes.get(idx1).add(e1);
        } else {
            System.out.println("node " + n1.name + " not found in graph");
        }
    }

    //----------------------------------------------------------------

    private void addResidualEdge(Node n1, Node n2, int capacity, boolean isBackward) {
        Edge e1 = new Edge(n1, n2, capacity, isBackward);
        int idx1 = this.nodes.indexOf(n1);
        if (idx1 >= 0) {
            this.nodes.get(idx1).addResidualEdge(e1);
        } else {
            System.out.println("node " + n1.name + " not found in graph");
        }
    }

    //----------------------------------------------------------------

    public void print() {
        for (Node n: this.nodes) {
            System.out.print("Node " + n.name + ":");
            for (Edge edge: n.adjlist) {
                System.out.print(" " + edge.n2.name + " (c=" + edge.capacity);
                System.out.print(", f=" + edge.flow + ")");
            }
            System.out.println();
        }
    }

    //----------------------------------------------------------------

    private void printResidual() {
        for (Node n: this.nodes) {
            System.out.print("Node " + n.name + ":");
            for (Edge edge: n.adjlistResid) {
                System.out.print(" " + edge.n2.name + " (c=" + edge.capacity);
                if (edge.isBackward)
                    System.out.print(" <=");
                System.out.print(")");
            }
            System.out.println();
        }
    }

    //----------------------------------------------------------------

    private List<Edge> findPathInResid(Node s, Node t) {
        int i, k, idx;
        boolean done, found;
        Node n1, n2;

        List<Edge> path = new ArrayList<Edge>();

        Stack<Node> stack = new Stack<Node>();
        boolean explored[] = new boolean[1 + this.nodes.size()];
        int parent[] = new int[1+this.nodes.size()];

        for (i=0; i<=this.nodes.size(); ++i)
            explored[i] = false;

        done = false;
        stack.push(s);
        while ( ! done && ! stack.empty() ) {
            n1 = stack.pop();
            if ( ! explored[n1.name] ) {
                explored[n1.name] = true;
                if (parent[n1.name] != 0)
                    System.out.println("tree: " + n1.name + " -> " + parent[n1.name]);
                for (Edge edge: n1.adjlistResid) {
                    n2 = edge.n2;
                    if ( ! explored[n2.name] ) {
                        stack.push(n2);
                        parent[n2.name] = n1.name;
                        if (n2.name == t.name)
                            done = true;
                    }
                }
            }
        }

        System.out.println("here's the backward path from " + t.name);
        done = false;
        idx = t.name;
        while ( ! done ) {
            if (parent[idx] == 0)
                done = true;
            else {
                System.out.println(parent[idx] + " to " + idx);
                // find the edge from parent[idx] to idx
                found = false;
                k = 0;
                while ( ! found && k < nodes.size()) {
                    if (nodes.get(k).name == parent[idx])
                        found = true;
                    else
                        k = k + 1;
                }
                n1 = nodes.get(k);
                found = false;
                for (k=0; ! found && k<n1.adjlistResid.size(); ++k) {
                    Edge e = n1.adjlistResid.get(k);
                    if (e.n2.name == idx) {
                        path.add(e);
                        found = true;
                    }
                }
                idx = parent[idx];
            }
        }

        System.out.println();
        return path;
    } // findPathInResid()

    //----------------------------------------------------------------

    public boolean checkFlow(Node s, Node t) {
        // check that flow out of s == flow into t
        // check conservation condition at each node
        boolean conservationCon = true; // assumed true until proven wrong
        boolean equal = true;
        boolean pass = true;
        boolean capacityCon = true;

        // print flow out of source
        // print flow into sink
        // print "checkFlow(): flow is : (valid or not)
        int sflow = 0;
        int tflow = 0;
        // check the flow leaving s and entering t
        for(int i = 0; i < s.adjlist.size(); i++) { // check flows leaving S
            sflow += s.adjlist.get(i).flow;
        }
        for(int i = 0; i < nodes.size(); i++){
            Node curr = nodes.get(i);
            for(int j = 0; j < curr.adjlist.size(); j++) {
                if (curr.adjlist.get(j).n2.name == t.name) {
                    tflow += curr.adjlist.get(j).flow;
                }
            }
        }
        System.out.print("flow out of source is: ");
        System.out.println(sflow);
        System.out.print("flow in to sink is: ");
        System.out.println(tflow);
        if(sflow != tflow){
            conservationCon = false;
        }

        List<Integer> flowIn = new ArrayList<Integer>(nodes.size() + 1);
        List<Integer> flowOut = new ArrayList<Integer>(nodes.size() + 1);
        for(int i = 0; i < flowIn.size(); i++){
            flowIn.set(i, 0);
            flowOut.set(i, 0);
        }

        // check the conservation conditions of the flow in the graph
        Node check;

        List<Edge> edges = getEdges();
        Edge curr;
        for(int i = 0; i < edges.size(); i++){
            curr = edges.get(i);
            flowIn.set(curr.n2.name,(flowIn.get(curr.n2.name) + curr.flow));
            flowOut.set(curr.n1.name, (flowOut.get(curr.n1.name) + curr.flow));
        }
        // remove first two elements and last since they are not being checked here
        flowIn.remove(0);
        flowIn.remove(1);
        flowIn.remove(flowIn.size()-1);
        flowOut.remove(0);
        flowOut.remove(1);
        flowOut.remove(flowIn.size()-1);

        for(int i = 0; i < flowIn.size(); i++){
            if(flowIn.get(i) != flowOut.get(i)){
                equal = false;
                System.out.println("checkFlow(): fail ");
                System.out.print("Node: ");
                System.out.print(flowIn.get(i));
                System.out.print("(flow out, flow in) : (");
                System.out.print(flowOut.get(i));
                System.out.print(", ");
                System.out.print(flowIn.get(i));
                System.out.println(")");
            }
        }
            // check for all nodes v in graph G that are not s and t that
            // flow entering the node is equal to the flow leaving the node
        if(equal && capacityCon && conservationCon){
            System.out.println("checkFlow(): flow is valid");
        } else {
            System.out.println("checkFlow(): flow is invalid");
            pass = false;
        }

        return pass; // only passable if all conditions are met
        // implement this
    }

    //----------------------------------------------------------------

    private void constructResidualGraph() {
        // implement this
        // create backwards edges of the current flow of the graph
        for(int i = 0; i < nodes.size(); i++){
            Node curr = nodes.get(i);
            curr.adjlistResid.clear();
        } // reset the previous residual graph to be recreated

        for(int i = 0; i < nodes.size(); i++){
            Node curr = nodes.get(i);
            for (int j = 0; j < curr.adjlist.size(); j++){
                // check flow and remaining capacity
                Edge currEdge = curr.adjlist.get(j);
                // check flow and create backward node from u, v in residual graph with
                if(currEdge.flow == 0) { // all free capacity
                    // create forward edge in residual for u to v with capacity and flow set to zero
                    Edge resid = new Edge(curr, currEdge.n2, currEdge.capacity, false);
                    curr.adjlistResid.add(resid); // add forward edge with the free capacity of the graph
                    System.out.println(currEdge.n1.name + " -> " + currEdge.n2.name + "(c=" + currEdge.capacity + ", f=" + currEdge.flow + "): f=0; create forward edge");
                } else { // check if capacity == flow
                    if (currEdge.flow == currEdge.capacity) { // if equal only a backwards edge is made
                        // backwards edge from u to v of the flow
                        Edge residadd = new Edge(currEdge.n2, currEdge.n1, currEdge.capacity, currEdge.flow, true);
                        currEdge.n2.adjlistResid.add(residadd);
                        System.out.println(currEdge.n1.name + " -> " + currEdge.n2.name + "(c=" + currEdge.capacity + ", f=" + currEdge.flow + "): f=c; create backward edge");
                    } else { // if not equal forward and backward are created
                        // forwards edge of the remaining capacity of the edge
                        Edge forward = new Edge(curr, currEdge.n2, (currEdge.capacity-currEdge.flow), false);
                        // backwards edge of the flow being transfered
                        Edge backward = new Edge(currEdge.n2, currEdge.n1, (currEdge.capacity - currEdge.flow), currEdge.flow, true);
                        currEdge.n2.adjlistResid.add(backward);
                        curr.adjlistResid.add(forward);
                        System.out.println(currEdge.n1.name + " -> " + currEdge.n2.name + "(c=" + currEdge.capacity + ", f=" + currEdge.flow + "): f<c; create forward edge");
                        System.out.println(currEdge.n1.name + " -> " + currEdge.n2.name + "(c=" + currEdge.capacity + ", f=" + currEdge.flow + "): f<c; create backward edge");

                    }
                }
            } // complete for all edges of node u
        } // complete for all nodes in graph G
    } // constructResidualGraph()

    //----------------------------------------------------------------

    public List<Edge> getEdges(){
        // return a list of the edges in the graph
        List<Edge> edges = new ArrayList<>();
        for(int i = 0; i < nodes.size(); i++){
            Node curr = nodes.get(i);
            for(int j = 0; j < curr.adjlist.size(); j++){
                edges.add(curr.adjlist.get(j));
            }
        }
        return edges;
    }

    //----------------------------------------------------------------

    private int findBottleneck(List<Edge> path) {
        // implement this
        // the bottle neck is the minimum residual capacity of any edge in P, with respect to the flow f
        // loop all residual edges for the minimum capacity
        int b = Integer.MAX_VALUE;
        for(int i = 0; i < path.size(); i++){
            Edge curr = path.get(i);
            if (curr.capacity < b) {
                b = curr.capacity;
            }
        }
        return b; // bottle neck value
    }

    //----------------------------------------------------------------

    private void augment(List<Edge> path) { // changes the flow of the graph related to te bottle neck
        // implement this
        // Let b = bottleneck(P, f) // f being the graph
        int b = findBottleneck(path);
        System.out.println("here's the path from s to t in Gf: ");
        for(int i = 0; i < path.size(); i++){
            System.out.println(path.get(i).n1.name + " to " + path.get(i).n2.name);
        }

        System.out.println("bottleneck = " + b);
        // for each edge (u, v) that is an elemet of P
        for (int i = 0; i < path.size(); i++) {
            // if e = (u, v) is a forward edge then
            Edge curr = path.get(i);
            if(curr.isBackward == false) { // fucked up
                // increase f(e) in G by b
                // only adjusting the value in the residual now fix for the acutal edge the residual represents
                Node adj = nodes.get(path.get(i).n1.name - 1); // minus 1 to get the node as it starts form 0
                // check that it adding to the edge that ends in n2 = path.n2
                for(int j = 0; j < adj.adjlist.size(); j++) {
                    if(adj.adjlist.get(j).n2 == path.get(i).n2){
                        adj.adjlist.get(j).flow += b;
                        System.out.println("forward edge " + path.get(i).n1.name + " -> " + path.get(i).n2.name + ": increase by " + b);
                    }
                }
            } else { // else ((u,v) is a backward edge, and let e = (v,u))
                // decrease f(e) in G by b
                Node adj = nodes.get(path.get(i).n1.name - 1);

                for(int j = 0; j < adj.adjlist.size(); j++) {
                    if(adj.adjlist.get(j).n2 == path.get(i).n2){
                        adj.adjlist.get(j).flow -= b;
                        System.out.println("backward edge " + path.get(i).n1.name + " -> " + path.get(i).n2.name + ": decreased by " + b);
                    }
                }
            }
            // for all edges on the path adjust the flow
        }
        System.out.println("here's the graph after augmenting");
        print();
    }

    //----------------------------------------------------------------

    public int maxFlow(Node s, Node t) {
        // implement this
        int flow = 0; // flow
        // while there is an s - t path in the residual graph (Gf)
        constructResidualGraph();
        printResidual();
        List<Edge> path = findPathInResid(s, t);
        while(path.size() > 0) {
            // let P be the simple path from s - t in Gf
            // f' = augment(f, P)
            augment(path);
            // update f to f'
            // update the residual graph Gf to be Gf'
            checkFlow(s, t);
            constructResidualGraph();
            printResidual();
            // refind path to see if a better solution is possible
            path = findPathInResid(s, t);

        }
        for(int i = 0; i < t.adjlistResid.size(); i++){
                flow += t.adjlistResid.get(i).flow;
        }

        System.out.println("max flow is " + flow);

        return flow;
    } // maxFlow()
}



