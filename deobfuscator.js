// 필요한 라이브러리 임포트
const fs = require('fs');
const beautify = require('js-beautify').js;
const esprima = require('esprima');

// 단계 1: 코드 포매팅 (Beautifying)
function beautifyCode(obfuscatedCode) {
  const beautifiedCode = beautify(obfuscatedCode, { indent_size: 2 });
  return beautifiedCode;
}

// 단계 2: 제어 흐름 분석 (Eval 사용 여부 체크)
function detectEvalUsage(code) {
  const evalRegex = /\beval\(/g;
  const match = code.match(evalRegex);
  if (match) {
    console.log("eval 함수 사용됨!");
  } else {
    console.log("eval 함수 사용되지 않음.");
  }
}

// 단계 3: 동적 값 추적
function trackDynamicValues(code) {
  const codeWithLogging = code.replace(/(var|let|const)\s+(\w+)/g, (match, keyword, varName) => {
    return `${match} // ${varName} 초기화\nconsole.log("${varName}:", ${varName});`;
  });
  console.log("변수 추적 코드:", codeWithLogging);
  return codeWithLogging;
}

// 단계 4: 코드 파싱 (AST 분석)
function parseCodeToAST(code) {
  try {
    const ast = esprima.parseScript(code);
    console.log("AST 구조:", JSON.stringify(ast, null, 2));
    return ast;
  } catch (e) {
    console.error("구문 분석 실패:", e);
  }
}

// 파일 읽기 (obfuscated_code.js)
function readObfuscatedCode(filePath) {
  return fs.readFileSync(filePath, 'utf-8');
}

// 파일 쓰기 (deobfuscated_code.js)
function writeDeobfuscatedCode(filePath, code) {
  fs.writeFileSync(filePath, code, 'utf-8');
}

// 메인 함수: 난독화된 파일을 읽고, 처리 후 복원된 코드를 새로운 파일에 저장
function processObfuscatedCode() {
  const obfuscatedFilePath = 'obfuscated_code.js';  // 입력 파일 경로
  const deobfuscatedFilePath = 'deobfuscated_code.js';  // 출력 파일 경로

  // 1단계: 난독화된 코드 읽기
  const obfuscatedCode = readObfuscatedCode(obfuscatedFilePath);
  console.log("난독화된 코드 읽기 완료");

  // 2단계: 코드 포매팅 (Beautifying)
  const beautifiedCode = beautifyCode(obfuscatedCode);
  console.log("코드 포매팅 완료");

  // 3단계: 제어 흐름 분석 (Eval 사용 여부 체크)
  detectEvalUsage(beautifiedCode);

  // 4단계: 동적 값 추적
  trackDynamicValues(beautifiedCode);

  // 5단계: 코드 파싱 (AST 분석)
  parseCodeToAST(beautifiedCode);

  // 6단계: 복원된 코드 deobfuscated_code.js로 저장
  writeDeobfuscatedCode(deobfuscatedFilePath, beautifiedCode);
  console.log(`${deobfuscatedFilePath}로 복원된 코드 저장 완료`);
}

// 메인 실행
processObfuscatedCode();
