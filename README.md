# JS-Deobfuscator

This project is a script that reads obfuscated JavaScript code, analyzes it, and restores it. The script formats the obfuscated code, checks for the usage of the `eval` function, tracks dynamic values, and analyzes the code through Abstract Syntax Tree (AST). Finally, the restored code is saved to a new file.

## Features

This script provides the following main features:

1. **Code Beautifying**: It reads obfuscated code and formats it in a readable way, making it easier to understand.
2. **Control Flow Analysis (Eval usage check)**: It checks whether the `eval()` function is used in the code.
3. **Dynamic Value Tracking**: It adds log statements to track variables during initialization, helping to identify dynamically changing values.
4. **Code Parsing (AST Analysis)**: It parses the code into an Abstract Syntax Tree (AST) and outputs the structure of the code for analysis.

## Requirements

To run this script, you will need Node.js and a few libraries.

### Required Libraries

- **js-beautify**: A library for formatting the code.
- **esprima**: A library for parsing JavaScript code.
- **fs**: A built-in library for working with the file system.

These libraries can be installed via `npm install`.

## Installation and Execution

### 1. Clone the Project

First, clone the project or download the files.

```bash
git clone https://github.com/korkeep/JS-deobfuscator.git
```

### 2. Install Required Libraries

Install the necessary libraries using the following command:

```bash
npm install js-beautify esprima
```

### 3. Prepare the Obfuscated Code File

Place the obfuscated JavaScript file in the project folder with the name `obfuscated_code.js`.

### 4. Run the Script

Execute the script with the following command to analyze the obfuscated code and generate the restored code:

```bash
node deobfuscator.js
```

Once the script runs, the obfuscated code will be formatted and analyzed, and the restored code will be saved in a new file named `deobfuscated_code.js`.

## Code Explanation

1. **Code Beautifying**: The `beautifyCode()` function formats the obfuscated code into a more readable format.
2. **Control Flow Analysis**: The `detectEvalUsage()` function checks for the usage of the `eval()` function in the code.
3. **Dynamic Value Tracking**: The `trackDynamicValues()` function adds log statements to track variables' values during initialization.
4. **Code Parsing**: The `parseCodeToAST()` function parses the JavaScript code into an Abstract Syntax Tree (AST) for analysis.

## File Description

- **obfuscated_code.js**: The file containing the obfuscated JavaScript code to be analyzed.
- **deobfuscated_code.js**: The file where the restored code will be saved.