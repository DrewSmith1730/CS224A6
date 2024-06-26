// jdh CS 224 Spring 2023

import java.util.List;

public class Main {
    public static void main(String args[]) {
        testOne();
        testTwo();
    }

    //---------------------------------------------------------------------
    // the example on p. 339, with s as node 1, u as node 2, v as node 3,
    // and t as node 4; and with initial flow set to zero

    public static void testOne() {
        Node n1 = new Node(1);
        Node n2 = new Node(2);
        Node n3 = new Node(3);
        Node n4 = new Node(4);

        Graph G = new Graph();
        G.addNode(n1);
        G.addNode(n2);
        G.addNode(n3);
        G.addNode(n4);

        // with flow
        G.addEdge(n1, n2, 20, 0);
        G.addEdge(n1, n3, 10, 0);
        G.addEdge(n2, n3, 30, 0);
        G.addEdge(n2, n4, 10, 0);
        G.addEdge(n3, n4, 20, 0);

        G.print();
        G.checkFlow(n1, n4);

        int mf = G.maxFlow(n1, n4);
        System.out.println("flow is " + mf);
        List<Edge> edges = G.getEdges();
        for (Edge e: edges) {
            System.out.println(e);
        }
    } // testOne()

    //--------------------------------------------------------------
    // example from the lecture notes, with initial flow of zero

    public static void testTwo() {
        Node n1 = new Node(1);
        Node n2 = new Node(2);
        Node n3 = new Node(3);
        Node n4 = new Node(4);
        Node n5 = new Node(5);

        Graph G = new Graph();
        G.addNode(n1);
        G.addNode(n2);
        G.addNode(n3);
        G.addNode(n4);
        G.addNode(n5);

        G.addEdge(n1, n2, 10,  4);
        G.addEdge(n1, n4,  3,  0);
        G.addEdge(n2, n3,  8,  4);
        G.addEdge(n2, n4,  5,  0);
        G.addEdge(n3, n4,  4,  4);
        G.addEdge(n3, n5,  2,  0);
        G.addEdge(n4, n5, 10,  4);

        G.print();
        G.checkFlow(n1, n5);

        G.maxFlow(n1, n5);
        G.print();
        G.checkFlow(n1, n5);
    } // testTwo()
}

