#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#define N 3000000          // 부모 DNA 시퀀스 길이
#define M 10000            // short read 개수
#define READ_LENGTH 32     // short read 길이
#define PARTS 4            // 돌연변이 삽입을 위해 분배한 파트 수수
#define MUTATIONS 3        // 돌연변이 개수 (=mismatch 개수)
#define NUM_CANDIDATES 101 // 아이 후보 수(진짜 1명 + 가짜 100명)

// 1) encode:  "A","C","G","T"를 -> 2bit 코드로 인코딩하는 함수
// 'A' -> 00
// 'C' -> 01
// 'G' -> 10
// 'T' -> 11
uint8_t encode(char dna) {
    switch (dna) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 0xFF; // 다른 입력 들어올 시 에러 처리
    }
}


// 2) decode: 2bit 코드를 -> 다시 ""A","C","G","T"로 디코딩하는 함수
// 00 -> 'A'
// 01 -> 'C'
// 10 -> 'G'
// 11 -> 'T'
char decode(uint8_t code) {
    switch (code) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N'; // 다른 입력 들어올 시 에러 처리
    }
}

// 3) generate_parentDNA: 부모 DNA 시퀀스 생성 함수
// : 랜덤
uint8_t* generate_parentDNA(int n) {
    uint8_t* seq = (uint8_t*)malloc(n); // dna sequence 저장
    const char list[] = "ATCG"; //염기 리스트
    int i;
    for (i = 0; i < n; i++) {
        char dna = list[rand() % 4];  // A,C,G,T 중 랜덤 선택
        seq[i] = encode(dna); // A,C,G,T -> 2bit로 인코딩딩해 저장장
    }
    return seq; // 부모 (전체) sequence 반환
}


// 4) generate_mutated :부모 DNA로부터 돌연변이가 삽입된 short read 생성 함수
// 부모 DNA sequence에서 시작위치를 기준으로 길이 32만큼 잘라
// 32 길이만큼 자른 부분을 4파트로 또 나눠 
// 4파트 중 3 파트에 랜덤하게 돌연변이 삽입해 아이 read 만들기
// -> 최종적으로 2bit 인코딩된 64bit read로 리턴턴
uint64_t generate_mutated_read(const uint8_t* parent_seq, int start_pos) {
    uint64_t read = 0;
    int i, part_size = READ_LENGTH / PARTS; // 32 / 4 = 8
    int part_idx, position, global_pos;
    uint8_t original, new_dna;
    uint64_t mask;

    // 원본 read 생성
    for (i = 0; i < READ_LENGTH; i++) {
        read <<= 2; // 2bit 왼쪽 이동(shift) (기존 염기 자리 이동)
        read |= parent_seq[start_pos + i]; // 2bit로 인코딩된 염기를 read에 추가
    }

    // MUTATIONS 만큼 돌연변이 생성
    for (i = 0; i < MUTATIONS; i++) {
        part_idx = rand() % PARTS; // 0 1 2 3 중 랜덤하게 파트 선택 후, 
        position = rand() % part_size; // 파트 안에서의 위치 선택
        global_pos = part_idx * part_size + position; //read 전체 기준 global 위치

        // 기존존 read에서 특정 위치(global_pos)의 염기를 가져오기 위해
        // 전체 read를 오른쪽으로 쉬프트하여 해당 위치가 가장 오른쪽에 오도록 이동
        // 이후 & 0x3을 수행해 해당 염기의 2bit만 추출
        original = (read >> (2 * (READ_LENGTH - 1 - global_pos))) & 0x3;
        // 새로운 염기서열을 랜덤하게 선택하되
        // 기존 염기서열과 중복되지는 않도록
        do {
            new_dna = rand() % 4;
        } while (new_dna == original);
        mask = ~(0x3ULL << (2 * (READ_LENGTH - 1 - global_pos)));
        read = (read & mask) | ((uint64_t)new_dna << (2 * (READ_LENGTH - 1 - global_pos)));
    }

    return read;
}


// 5) generate_fake_read :완전 랜덤 fake read 생성 함수
// 아이 후보지만, 정작 친자식이 아닌 경우
uint64_t generate_fake_read() {
    uint64_t read = 0;
    int i;
    for (i = 0; i < READ_LENGTH; i++) {
        read <<= 2;  // 왼쪽으로 2bit 이동시켜 자리 비우기
        read |= (rand() % 4);  // 0~3(염기 A, C, G, T)을 랜덤하게 선택 후 삽입
        // rand() % 4 값은 넷 중 하나
        // 0 : A
        // 1 : C
        // 2 : G
        // 3 : T

    }
    return read; // 최종적으로 2bit 단위로 32개의 염기가 저장된 read 반환
}


// txt함수) save_parentDNA : 부모 DNA를 parentDNA.txt 파일로 저장하는 용도도
void save_parentDNA(const uint8_t* seq, int length) {
    FILE* fp = fopen("parentDNA.txt", "w");
    if (!fp) {
        perror("파일 열기 실패 (parentDNA.txt)");
        return;
    }
    int i;
    for (i = 0; i < length; i++) {
        fprintf(fp, "%c", decode(seq[i]));
    }
    fclose(fp);
}


// txt함수) txt_ childReads : child_shortReads.txt 저장 함수
void save_child_reads(uint64_t* reads, int m) {
    FILE* fp = fopen("child_shortReads.txt", "w");
    if (!fp) {
        perror("파일 열기 실패 (child_shortReads.txt)");
        return;
    }
    int i, j;
    uint8_t dna;
    for (i = 0; i < m; i++) { //각 read마다 반복
        for (j = 0; j < READ_LENGTH; j++) { 
            // 64비트 read에서 상위 비트부터 차례로 2비트씩 잘라 디코딩
            dna = (reads[i] >> (2 * (READ_LENGTH - 1 - j))) & 0x3;
            fprintf(fp, "%c", decode(dna));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}


// txt함수) childDNA.txt 저장 함수
void save_child_seq(const uint8_t* parent_seq, int length) {
    FILE* fp = fopen("childDNA.txt", "w");
    if (!fp) {
        perror("파일 열기 실패 (childDNA.txt)");
        return;
    }
    int i;
    for (i = 0; i < length; i++) {
        fprintf(fp, "%c", decode(parent_seq[i]));
    }
    fclose(fp);
}


// 6) reconstructing: short reads 정보를 기반으로 전체 DNA 시퀀스 재구성
// 입력 (매개변수)
//      reads: short read 배열(각 read는 64bit로 표현)
//      positions : 각 read의 시작 위치 배열
//      m : 한 아이로부터 추출되는 short read 개수
// 출력 (반환값)
//      reconstruct된 시퀀스 배열((uint8_t*), 각 염기가 2bit인코딩
uint8_t* reconstructing(uint64_t* reads, int* positions, int m) {
    uint8_t* seq = (uint8_t*)malloc(N); // 최종
    memset(seq, 0xFF, N); // 0xFF: 'N'으로 초기값 설정 

    int* covered = (int*)calloc(N, sizeof(int)); // 각 위치에 얼마나 read가 덮였는지 세기
    int** dna_counts = (int**)malloc(N * sizeof(int*)); // 각 위치별 ACGT 카운트 배열
    int i, j, k, start_pos, idx, max;

    // 4개(A,C,G,T) 카운팅 배열 초기화
    for (i = 0; i < N; i++) {
        dna_counts[i] = (int*)calloc(4, sizeof(int)); //각 위치에 4개 염기 카운터 생성
    }

    // 각 read의 dna 카운트
    for (i = 0; i < m; i++) { // 모든 read 반복
        start_pos = positions[i];  // read의 시작위치를 정하고
        for (j = 0; j < READ_LENGTH; j++) { // read의 각 염기마다 반복하며
            if (start_pos + j >= N) continue; //boundary 체크
            uint8_t dna = (reads[i] >> (2 * (READ_LENGTH - 1 - j))) & 0x3; // 2bit씩 이동(shift)하며 mask로 염기 추출
            dna_counts[start_pos + j][dna]++; // 해당 위치에 해당 염기 카운트 증가
            covered[start_pos + j]++;  // covered 증가
        }
    }

    // majority vote을 활용해 reconstruct
    for (i = 0; i < N; i++) {
        if (covered[i] == 0) {
            seq[i] = 0xFF;   // read로 덮이지 않았으면 'N'으로 표시
        } else {
            max = 0;
            idx = 0;
            for (k = 0; k < 4; k++) { // A,C,G,T 중 가장 많은 염기를 선택해
                if (dna_counts[i][k] > max) {
                    max = dna_counts[i][k];
                    idx = k;
                }
            }
            seq[i] = idx; // 선택된 염기를 reconstruct된 시퀀스에 삽입
        }
    }

    // (동적) 메모리 해제하기
    for (i = 0; i < N; i++) {
        free(dna_counts[i]);
    }
    free(dna_counts);
    free(covered);

    return seq;
}


// 7) calc_similarity: parent_seq와 child_seq의 유사도를 계산 함수
// : 부모 DNA 시퀀스와 majority vote(최빈값)를 통해 reconstruct한 아이의 dna sequence 유사도 비교하기
// 수치 -> %로 표현
double calc_similarity(const uint8_t* parent_seq, const uint8_t* child_seq) {
    int i, match = 0, total = 0;
    for (i = 0; i < N; i++) {
        if (child_seq[i] == 0xFF) continue;    // child_seq가 비어있으면 건너뛰기
        if (parent_seq[i] == child_seq[i]) match++; // 동일하면 match++
        total++;
    }

    double similarity;
    if (total == 0) {
        similarity = 0.0;
    } else {
        similarity = ((double)match / total) * 100.0;
    }

    return similarity; //유사도 반환
}

// main 함수------------------------------------------------------------------------------
int main() {
    int i, j, pos, step;
    srand(time(NULL)); //난수 초기화

    printf("--- 데이터 소스 생성 시작 ---\n");

    // 1. 부모 DNA 시퀀스 생성하고 .txt로 저장하기
    uint8_t* parent_seq = generate_parentDNA(N);
    printf("부모 시퀀스 생성 완료.\n\n");
    save_parentDNA(parent_seq, N);
    printf("parentDNA.txt 저장 완료.\n");

    // 2. 후보 배열 선언하기 (candidates array)
    uint64_t* cand_reads[NUM_CANDIDATES]; //short reads 저장하는 배열
    int* cand_pos[NUM_CANDIDATES]; // 각 read의 시작 위치
    char* cand_ids[NUM_CANDIDATES]; // 후보 이름 (id)
    double cand_similarity[NUM_CANDIDATES]; // 후보별 유사도
    uint8_t* cand_DNAseq[NUM_CANDIDATES]; // majority vote방식으로 복원된 seq

    // 3. 친자(==아이후보 1) short reads생성하고 .txt파일로 저장
    cand_ids[0] = strdup("아이 후보 1 "); //후보 이름(id)
    cand_reads[0] = (uint64_t*)malloc(sizeof(uint64_t) * M);
    cand_pos[0] = (int*)malloc(sizeof(int) * M);
    step = (N - READ_LENGTH) / M; //read 간격 계산하기
    for (i = 0; i < M; i++) {
        pos = i * step;
        cand_reads[0][i] = generate_mutated_read(parent_seq, pos);//변이를 집어넣어 short read 생성하기
        cand_pos[0][i] = pos;
    }
    save_child_reads(cand_reads[0], M); // short read .txt로 저장
    printf("child_shortReads.txt 저장 완료.\n");

    // 4. child_seq .txt로 저장하기
    uint8_t* child_seq = (uint8_t*)malloc(N);
    memcpy(child_seq, parent_seq, N); //memcpy로 부모 DNA 복사하기
    save_child_seq(child_seq, N);     //childDNA.txt로 저장
    printf("childDNA.txt 저장 완료.\n");

    // 5. (친자가 아닌 fake)아이 후보  2~101명 정도의 short reads 생성
    printf("\n--- 아이 후보 분석 시작 ---\n");
    for (i = 1; i < NUM_CANDIDATES; i++) {
        char label[50];
        sprintf(label, "아이 후보 %d ", i + 1); //label 재열에 "아이 후보 %d(id)"를 문자열로 만들어 넣기 
        cand_ids[i] = strdup(label); // label에 들어있는 문자열 "아이 후보 id"를 힙에 동적으로 복사 ( 새로운 문자열 메모리 만들어 리턴)
                                     // *label이 함수 안에서 사라져도 새로운 문자열 메모리는 계속 유지되도록 
        cand_reads[i] = (uint64_t*)malloc(sizeof(uint64_t) * M);
        cand_pos[i] = (int*)malloc(sizeof(int) * M);
        for (j = 0; j < M; j++) {
            cand_reads[i][j] = generate_fake_read(); //랜덤 fake read만들기
            cand_pos[i][j] = rand() % (N - READ_LENGTH + 1); //랜덤 위치
        }
    }

    // 6. 결과(result) 출력 및 .txt파일에 기록하기
    FILE* result_fp = fopen("result.txt", "w");
    if (!result_fp) {
        perror("파일 열기 실패 (result.txt)");
        return 1;
    }

    // 후보별 유사도 분석하기 & 출력
    for (i = 0; i < NUM_CANDIDATES; i++) {
        cand_DNAseq[i] = reconstructing(cand_reads[i], cand_pos[i], M);
        cand_similarity[i] = calc_similarity(parent_seq, cand_DNAseq[i]);
        printf("%s 최종 유사도: %.4f%%\n", cand_ids[i], cand_similarity[i]); //소수점까지 파악하기 위해 .4f% 사용
        fprintf(result_fp, "%s 최종 유사도: %.4f%%\n", cand_ids[i], cand_similarity[i]);
    }

    // 7. 최고 유사도 후보 판별
    int best_idx = 0;
    double most_similar = cand_similarity[0];
    for (i = 1; i < NUM_CANDIDATES; i++) {
        if (cand_similarity[i] > most_similar) {
            most_similar = cand_similarity[i];
            best_idx = i;
        }
    }

    // 8. 최종 판별 결과 출력
    fprintf(result_fp, "\n=== 최종 판별 결과 ===\n");
    if (best_idx == 0) {
        printf("=> 최종 판별: 아이 %d명 중 아이 후보 %d번이 %.4f%%의 가장 높은 유사도로 친자식일 가능성이 가장 높습니다.\n",
               NUM_CANDIDATES, best_idx + 1, most_similar);
        fprintf(result_fp, "=> 최종 판별: 아이 %d명 중 아이 후보 %d번이 %.4f%%의 가장 높은 유사도로 친자식일 가능성이 가장 높습니다.\n",
                NUM_CANDIDATES, best_idx + 1, most_similar);
    } else {
        printf("=> 최종 판별: 이번 분석에서는 친자식을 찾지 못했습니다.\n");
        fprintf(result_fp, "=> 최종 판별: 이번 분석에서는 친자식을 찾지 못했습니다.\n");
    }

    fclose(result_fp);

    // 9. 메모리 해제
    free(parent_seq);
    free(child_seq);
    for (i = 0; i < NUM_CANDIDATES; i++) {
        free(cand_reads[i]);
        free(cand_pos[i]);
        free(cand_ids[i]);
        free(cand_DNAseq[i]);
    }

    printf("result.txt 저장 완료.\n");
    return 0;
}
