#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
using namespace std;

struct SequenceNode {
    string dnaSequence;
    int mutationCount;
    int sequenceIndex;
    SequenceNode* left;
    SequenceNode* right;

    SequenceNode(string seq, int mut, int index) {
        dnaSequence = seq;
        mutationCount = mut;
        sequenceIndex = index;
        left = nullptr;
        right = nullptr;
    }
};

class DNAAnalyzer {
private:
    string referenceSequence;
    SequenceNode* root;

    string convertToUpperCase(string seq) {
        for (int i = 0; i < seq.length(); i++) {
            seq[i] = toupper(seq[i]);
        }
        return seq;
    }

    bool isValidDNASequence(string seq) {
        for (int i = 0; i < seq.length(); i++) {
            if (seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'G' && seq[i] != 'C' && seq[i] != '-') {
                return false;
            }
        }
        return true;
    }

    SequenceNode* insertInBST(SequenceNode* node, string seq, int mutations, int index) {
        if (node == nullptr) {
            return new SequenceNode(seq, mutations, index);
        }

        if (mutations < node->mutationCount) {
            node->left = insertInBST(node->left, seq, mutations, index);
        } else {
            node->right = insertInBST(node->right, seq, mutations, index);
        }

        return node;
    }

    void displayInOrder(SequenceNode* node, int& rank) {
        if (node == nullptr) {
            return;
        }

        displayInOrder(node->left, rank);

        cout << "Rank " << rank << ":" << endl;
        cout << "Sequence #" << node->sequenceIndex << ": " << node->dnaSequence << endl;
        cout << "Mutations: " << node->mutationCount << endl;
        cout << "Similarity: " << fixed << setprecision(2)
             << (1.0 - (double)node->mutationCount / referenceSequence.length()) * 100 << "%" << endl;
        cout << "--------------------------------" << endl;

        rank++;

        displayInOrder(node->right, rank);
    }

    void printTree(SequenceNode* node, int level = 0, string prefix = "Root: ") {
        if (node != nullptr) {
            cout << string(level * 4, ' ') << prefix << "Seq #" << node->sequenceIndex 
                 << " (Mutations: " << node->mutationCount << ")" << endl;
            printTree(node->left, level + 1, "L--- ");
            printTree(node->right, level + 1, "R--- ");
        }
    }

    void deleteBST(SequenceNode* node) {
        if (node == nullptr) {
            return;
        }

        deleteBST(node->left);
        deleteBST(node->right);
        delete node;
    }

public:
    DNAAnalyzer(string refSeq) {
        refSeq = convertToUpperCase(refSeq);

        while (!isValidDNASequence(refSeq)) {
            cout << "Error: Reference sequence must contain only A, T, G, C, and -. Please re-enter: ";
            cin >> refSeq;
            refSeq = convertToUpperCase(refSeq);
        }

        referenceSequence = refSeq;
        root = nullptr;
    }

    int countMutations(string seq) {
        if (seq.length() != referenceSequence.length()) {
            return -1;
        }

        int mutationCount = 0;
        for (int i = 0; i < seq.length(); i++) {
            if (toupper(seq[i]) != referenceSequence[i]) {
                mutationCount++;
            }
        }
        return mutationCount;
    }

    void addNewSequence(string seq, int index) {
        seq = convertToUpperCase(seq);

        while (!isValidDNASequence(seq) || countMutations(seq) == -1) {
            if (!isValidDNASequence(seq)) {
                cout << "Error: Sequence must contain only A, T, G, C and -. Please re-enter: ";
            } else {
                cout << "Error: Sequence must be of equal length (" << referenceSequence.length() << " characters). Please re-enter: ";
            }
            cin >> seq;
            seq = convertToUpperCase(seq);
        }

        int mutations = countMutations(seq);
        root = insertInBST(root, seq, mutations, index);
    }

    void showAnalysisResults() {
        cout << "\nReference sequence: " << referenceSequence << endl;
        
        cout << "\nDNA Mutation Tree Structure:" << endl;
        cout << "--------------------------------" << endl;
        printTree(root);
        cout << "--------------------------------" << endl;
        
        cout << "\nAnalysis Results (sorted by relevance):" << endl;
        cout << "--------------------------------" << endl;

        int rank = 1;
        displayInOrder(root, rank);
    }

    ~DNAAnalyzer() {
        deleteBST(root);
    }
};

struct AlzheimerMutationNode {
    string mutationSequence;
    string relatedGene;
    AlzheimerMutationNode *left, *right;

    AlzheimerMutationNode(string seq, string geneName) {
        mutationSequence = seq;
        relatedGene = geneName;
        left = right = nullptr;
    }
};

AlzheimerMutationNode* addMutationSequence(AlzheimerMutationNode* root, string seq, string geneName) {
    if (!root) return new AlzheimerMutationNode(seq, geneName);
    if (seq < root->mutationSequence)
        root->left = addMutationSequence(root->left, seq, geneName);
    else if (seq > root->mutationSequence)
        root->right = addMutationSequence(root->right, seq, geneName);
    return root;
}

void printAlzheimerTree(AlzheimerMutationNode* node, int level = 0, string prefix = "Root: ") {
    if (node != nullptr) {
        cout << string(level * 4, ' ') << prefix << node->mutationSequence 
             << " (" << node->relatedGene << ")" << endl;
        printAlzheimerTree(node->left, level + 1, "L--- ");
        printAlzheimerTree(node->right, level + 1, "R--- ");
    }
}

vector<string> findGenesInSequence(AlzheimerMutationNode* root, string seq) {
    vector<string> matchedGenes;
    if (!root) return matchedGenes;

    if (seq.find(root->mutationSequence) != string::npos) {
        matchedGenes.push_back(root->relatedGene);
    }

    vector<string> leftMatches = findGenesInSequence(root->left, seq);
    vector<string> rightMatches = findGenesInSequence(root->right, seq);
    
    matchedGenes.insert(matchedGenes.end(), leftMatches.begin(), leftMatches.end());
    matchedGenes.insert(matchedGenes.end(), rightMatches.begin(), rightMatches.end());

    return matchedGenes;
}

void detectAlzheimersGenes(string refSeq, vector<string> sequences) {
    AlzheimerMutationNode* mutationRoot = nullptr;

    mutationRoot = addMutationSequence(mutationRoot, "ATGCCGTGCGCG", "APOE");
    mutationRoot = addMutationSequence(mutationRoot, "ATTGAAAAATGC", "APP");
    mutationRoot = addMutationSequence(mutationRoot, "ATGCTCGGA", "PSEN1 (M146L)");
    mutationRoot = addMutationSequence(mutationRoot, "CTGGTCTTC", "PSEN1 (L286V)");
    mutationRoot = addMutationSequence(mutationRoot, "AACATTGGC", "PSEN2");

    cout << "\n--- ALZHEIMER'S SEQUENCE CHECK ---\n";
    
    cout << "\nAlzheimer's Gene Mutation Tree:" << endl;
    cout << "--------------------------------" << endl;
    printAlzheimerTree(mutationRoot);
    cout << "--------------------------------" << endl;

    string refUpper = refSeq;
    for (char& c : refUpper) c = toupper(c);

    vector<string> refGenes = findGenesInSequence(mutationRoot, refUpper);
    if (!refGenes.empty()) {
        cout << "\nREF SEQUENCE: INDICATES ALZHEIMER (MUTATIONS IN GENE ";
        for (string& gene : refGenes) cout << gene << " )";
        cout << endl;
    } else {
        cout << "\nREF SEQUENCE: IS NORMAL\n";
    }

    for (int i = 0; i < sequences.size(); i++) {
        string userSeq = sequences[i];
        for (char& c : userSeq) c = toupper(c);

        vector<string> matchedGenes = findGenesInSequence(mutationRoot, userSeq);
        if (!matchedGenes.empty()) {
            cout << "SEQ " << i + 1 << ": INDICATES ALZHEIMER MUTATIONS IN GENE(S): ";
            for (string& gene : matchedGenes) cout << gene << " ";
            cout << endl;
        } else {
            cout << "SEQ " << i + 1 << ": IS NORMAL\n";
        }
    }
}

int main() {
    string referenceInput;
    cout << "Enter reference DNA sequence: ";
    cin >> referenceInput;

    DNAAnalyzer analyzer(referenceInput);

    int totalSequences;
    cout << "Enter number of DNA sequences to compare: ";
    cin >> totalSequences;

    vector<string> userSequences;

    for (int i = 0; i < totalSequences; i++) {
        string seqInput;
        cout << "Enter sequence " << i + 1 << ": ";
        cin >> seqInput;
        userSequences.push_back(seqInput);
        analyzer.addNewSequence(seqInput, i + 1);
    }

    analyzer.showAnalysisResults();

    detectAlzheimersGenes(referenceInput, userSequences);

    return 0;
}