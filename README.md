# 🧬 Alzheimer’s Gene Mutation Detection using Binary Search Tree (BST)

## 📌 Overview

This project implements a **mutation-based Alzheimer’s detection system** using a **Binary Search Tree (BST)** for biologically ordered sequence comparison and ranking.

Unlike deep learning black-box models, this system provides:

* ✅ Interpretable mutation detection
* ✅ Gene-level traceability
* ✅ Similarity-based mutation ranking
* ✅ Real-time command-line feedback

The system focuses on detecting mutations in Alzheimer’s-associated genes using logical and structured sequence comparison.

---

# 🎯 Problem Motivation

Recent Alzheimer’s genomic studies rely heavily on deep learning models such as:

* CNN
* RNN
* Bi-LSTM with Attention
* Ensemble Machine Learning

Although they achieve high accuracy (up to 93%), they suffer from:

* ❌ High computational cost
* ❌ Poor interpretability
* ❌ Black-box behavior
* ❌ Lack of gene-level mutation traceability

This project addresses those limitations by implementing a **structured and interpretable tree-based mutation detection approach**.

---

# 🧠 Core Idea of This Implementation

Instead of training a neural network, this system:

1. Takes user DNA sequence input
2. Compares it against Alzheimer’s-associated genes
3. Detects base-level mutations
4. Stores mutation scores in a Binary Search Tree
5. Ranks sequences based on mutation similarity
6. Displays gene-level pattern matches

---

# 🏗️ System Architecture

```
User Input DNA Sequence
        ↓
Mutation Detection Engine
        ↓
Gene Pattern Matching
        ↓
BST-Based Ranking
        ↓
Similarity Score Output
```

---

# 🌳 Why Binary Search Tree (BST)?

The BST ensures:

* Biological order preservation
* Efficient mutation score insertion
* Sorted ranking of mutation similarity
* Fast retrieval of most similar gene match

This makes the system:

✔ Lightweight
✔ Fast
✔ Interpretable
✔ Structurally organized

---

# 📊 Comparison with Previous Research

| Approach                   | Accuracy      | Problem                  | Our Improvement                  |
| -------------------------- | ------------- | ------------------------ | -------------------------------- |
| CNN/RNN (2022)             | 92.45%        | Black-box, high cost     | Fully interpretable logic        |
| Bi-LSTM + Attention (2021) | 93.7%         | Poor gene-specific focus | Gene-level traceability          |
| Ensemble ML (2020)         | 88% Precision | Heavy preprocessing      | Direct mutation comparison       |
| Hash Tree (2023)           | 0.2s speed    | No biological order      | Ordered BST + similarity ranking |

---

# 🚀 Key Features

* 🧬 Gene-level Alzheimer’s mutation detection
* 🌳 BST-based mutation similarity ranking
* ⚡ Real-time CLI feedback
* 🔍 Base-level mutation comparison
* 📈 Ranked mutation similarity scoring
* 🧠 Fully interpretable logic

---

# 📂 Project Structure

```
Alzheimer-Mutation-Detection/
│
├── main.py
├── bst.py
├── mutation_detector.py
├── gene_database.txt
├── README.md
```

---

# 🛠 Technologies Used

* Python
* Data Structures (Binary Search Tree)
* Genomic Sequence Processing
* String Comparison Algorithms

---

# ▶️ How to Run

```bash
python main.py
```

Enter DNA sequence when prompted.

The system will:

* Detect mutation differences
* Compare against Alzheimer-associated genes
* Rank similarity
* Display results

---

# 🧪 Example Output

```
Detected Gene: APOE
Mutation Score: 7
Similarity Rank: 2nd Closest Match
Interpretation: Moderate similarity with known Alzheimer mutation pattern
```

---

# 📈 Advantages of This System

✔ No heavy computational training
✔ No synthetic encoding required
✔ Direct base-level mutation comparison
✔ Clear and explainable results
✔ Fast execution
✔ Structured tree-based ranking

---

# ⚠️ Current Limitations

* CLI-based (No GUI)
* No deep learning prediction probability
* Works on predefined Alzheimer gene patterns

---

# 🔮 Future Improvements

* Add graphical interface
* Integrate optional ML prediction module
* Expand gene database
* Add visualization of mutation sites
* Deploy as web-based tool

---

# 🎓 Academic Relevance

This project is suitable for:

* Data Structures & Algorithms coursework
* Bioinformatics research
* Mutation detection studies
* Explainable AI comparison studies
* Research proposal demonstrations

---

# 👨‍💻 Author

Muhammad Hammad
BS Bioinformatics
AI & Computational Biology Enthusiast

---

# 🌍 Future Vision

This project aims to combine:

* Data Structures
* Genomics
* Interpretable AI
* Research-oriented computing

to build explainable and efficient biomedical decision systems.

---

