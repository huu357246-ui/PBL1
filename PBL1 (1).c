#include <stdio.h>
#include <string.h>
#define eps 1e-12
#define MAX 100


typedef struct {
    double data[MAX];
    int n;
} Danhsachdac;

typedef struct {
    double A;
    double B;
    double a;
    double b;
} ketqua;


int readFile(char *namefile, Danhsachdac *dsx, Danhsachdac *dsy);
void out(Danhsachdac x, Danhsachdac y);
double Abs(double x);
double Pow(double x, int n);
double Epx(double x);
double lnx(double x);
double PowReal(double x, double n);
void saveFile(int n, ketqua kq, Danhsachdac x, Danhsachdac y);
int solve(double x[], double y[], int n, ketqua *kq);
void nhaptenfile(char *nameFile);
void check(int n, ketqua *kq, Danhsachdac x, Danhsachdac y);
void menuin(int luachon1, Danhsachdac x, Danhsachdac y, ketqua kq, int n);
void menu();


int main() {
    printf("----------------------------------------------------------------------------------\n");
    printf("| %-11s | %-64s |\n", "PBL1",        "DU AN LAP TRINH TINH TOAN");
    printf("----------------------------------------------------------------------------------\n");
    printf("| %-11s | %-64s |\n", "SINH VIEN",   "TRAN VAN ANH BINH");
    printf("| %-11s | %-64s |\n", "",             "NGUYEN TRONG HUU");
    printf("----------------------------------------------------------------------------------\n");
    printf("| %-11s | %-64s |\n", "DE TAI: 910", "TIM HE SO A, B CUA PHUONG TRINH Y=A*X^B");
    printf("| %-11s | %-64s |\n", "",             "BANG PHUONG PHAP BINH PHUONG BE NHAT");
    printf("----------------------------------------------------------------------------------\n");
    menu();
    return 0;
}


int readFile(char *namefile, Danhsachdac *dsx, Danhsachdac *dsy) {
    FILE *f = fopen(namefile, "r");
    if (f == NULL) {
        printf("LOI : Khong mo duoc file!");
        return 0;
    }
    int n;
    if (fscanf(f, "%d", &n) != 1) {
        printf("LOI: FILE khong chua gia tri n!");
        fclose(f);
        return 0;
    }
    if (n <= 0) {
        printf("LOI: n = %d khong hop le!\n", n);
        fclose(f);
        return 0;
    }

    if (n > MAX) {
        printf("LOI: n = %d vuot qua gia tri toi da la %d\n", n, MAX);
        fclose(f);
        return 0;
    }
    dsx->n = n;
    dsy->n = n;
    for (int i = 0; i < n; i++) {
        if (fscanf(f, "%lf %lf", &dsx->data[i], &dsy->data[i]) != 2) {
            printf("LOI: Doc du lieu tai dong %d that bai\n", i+1);
            fclose(f);
            return 0;
        }
    }
    fclose(f);
    return n;
}

void out(Danhsachdac x, Danhsachdac y) {
    printf("%-10s|%-15s|%-15s\n", "Diem", "x", "y");
    printf("----------+---------------+---------------\n");
    for (int i = 0; i < x.n; i++) {
        printf("%-10d|%15.6lf|%15.6lf\n", i + 1, x.data[i], y.data[i]);
    }
}

double Abs(double x) {
    if (x < 0) return -x;
    return x;
}

double Pow(double x, int n) {
    double p = 1;
    if (n >= 0) {
        for (int i = 0; i < n; i++) {
            p *= x;
        }
    }
    else if (n < 0) {
        for (int i = 0; i < -n; i++) {
            p /= x;
        }
    }
    return p;
}

double Epx(double x) {
    if (Abs(x)>500){
        return 0;
    }
    double s = 1.0;
    double a = 1.0;
    int i = 1;
    while (Abs(a) > eps) {
        a = a * (x / i);
        s += a;
        i++;
    }
    return s;
}

double lnx(double x) {
    double s = 0, m = (x - 1.0) / (x + 1.0), n = m * m, p = m;
    int i = 1;
    do {
        s += p / i;
        p *= n;
        i += 2;
    } while (Abs(p / i) > eps);
    return 2 * s;
}

double PowReal(double x, double n) {
    return Epx(n * lnx(x));
}

void saveFile(int n, ketqua kq, Danhsachdac x, Danhsachdac y) {
    char nameFile[256];
    printf("Nhap ten file de luu ket qua: ");
    scanf("%255s", nameFile);
    strncat(nameFile, ".txt", 255 - strlen(nameFile));

    FILE *f = fopen(nameFile, "w");
    if (f == NULL) {
        printf("LOI: Khong the tao file %s!\n", nameFile);
        return;
    }

    fprintf(f, "=== KET QUA TINH TOAN ===\n");
    fprintf(f, "Phuong trinh: y = %.6lf * x^%.6lf\n", kq.a, kq.b);
    fprintf(f, "  a = e^%.6lf = %.6lf\n", kq.A, kq.a);
    fprintf(f, "  b = %.6lf\n\n", kq.b);

    fclose(f);
    printf("Da luu ket qua vao file %s\n", nameFile);
}

int solve(double x[], double y[], int n, ketqua *kq) {
    double TongX = 0, TongY = 0, TongXY = 0, TongX2 = 0;
    printf("%-6s|%-12s|%-12s|%-12s|%-12s\n", "Buoc", "Tong X", "Tong Y", "Tong XY", "Tong X^2");
    printf("------+------------+------------+------------+------------\n");
    for (int i = 0; i < n; i++) {
        double lnX = lnx(x[i]);
        double lnY = lnx(y[i]);
        if (x[i] < 0) {
            printf("Gia tri x[%d] = %lf khong hop le, khong the tinh ln(x).\n", i + 1, x[i]);
            return 0;
        }
        if (y[i] < 0) {
            printf("Gia tri y[%d] = %lf khong hop le, khong the tinh ln(y).\n", i + 1, y[i]);
            return 0;
        }
        TongX  += lnX;
        TongY  += lnY;
        TongXY += lnX * lnY;
        TongX2 += Pow(lnX, 2);
        printf("%-6d|%12.6lf|%12.6lf|%12.6lf|%12.6lf\n", i + 1, TongX, TongY, TongXY, TongX2);
    }
    printf("=> Tong X = %lf\n", TongX);
    printf("   Tong Y = %lf\n", TongY);
    printf("   Tong XY = %lf\n", TongXY);
    printf("   Tong X^2 = %lf\n", TongX2);
    double mau = (n * TongX2 - Pow(TongX, 2));
    if (Abs(mau) < eps) return 0;
    printf("Ta co he phuong trinh:\n");
    printf("%d*A + %lf*B = %lf (1)\n", n, TongX, TongY);
    printf("%lf*A + %lf*B = %lf (2)\n", TongX, TongX2, TongXY);
    kq->B = (n * TongXY - TongX * TongY) / mau;
    kq->A = (TongY - kq->B * TongX) / n;
    printf("=> A = (%lf - %lf*%lf)/%lf = %lf\n", TongY, kq->B, TongX, n, kq->A);
    printf("   B = (%d*%lf - %lf*%lf)/(%d*%lf - %lf^2) = %lf\n", n, TongXY, TongY, TongX, n, TongX2, TongX, kq->B);
    if (Epx(kq->A) == 0) {
        printf("Gia tri cua a qua nho, khong the tinh duoc\n");
        return 0;
    }
    kq->a = Epx(kq->A);
    kq->b = kq->B;
    return 1;
}

void nhaptenfile(char *nameFile) {
    printf("Nhap ten file: ");
    scanf("%255s", nameFile);
    strncat(nameFile, ".txt", 255 - strlen(nameFile));
}

void check(int n, ketqua *kq, Danhsachdac x, Danhsachdac y) {
    double epscheck = 0.05;  // sai số tương đối 0.01%
    for (int i = 0; i < n; i++) {
        double ytinh = kq->a * PowReal(x.data[i], kq->b);
        double saisotuongdoi = Abs(y.data[i] - ytinh) / Abs(y.data[i]);
        if (saisotuongdoi > epscheck) {
            printf("Ket qua khong chinh xac tai diem (%lf, %lf)\n", x.data[i], y.data[i]);
        } else {
            printf("Ket qua chinh xac tai diem (%lf, %lf)\n", x.data[i], y.data[i]);
        }
    }
}

void menuin(int luachon1, Danhsachdac x, Danhsachdac y, ketqua kq, int n) {
    int luachon2, t_check=0;
    do {
        printf("===================================== MENU =====================================\n");
        printf("1. Xuat du lieu cac toa do (x, y).\n");
        printf("2. Tinh he so a, b.\n");
        printf("3. Kiem tra ket qua tai cac diem du lieu.\n");
        printf("4. Luu ket qua vao file.\n");
        printf("5. Quay lai.\n");
        printf("Lua chon: ");
        scanf("%d", &luachon2);
        switch (luachon2)
        {
        case 1:
            out(x, y);
            break;

        case 2:
            if(n<2){
                printf("Phai co it nhat hai diem du lieu!\n");
            } else if (n>=2){
                if (solve(x.data, y.data, n, &kq)) {
                    t_check=1;
                    printf("Gia tri cua a = e^%lf = %lf.\n", kq.A, kq.a);
                    printf("Gia tri cua b = %lf.\n", kq.b);
                    printf("Phuong trinh ban dau la: y = %lf * x^%lf.\n", kq.a, kq.b);
                } else {
                    printf("Khong the tinh duoc he so (mau = 0)!\n");
            }
        }
            break;
        case 3:
            if (t_check) {
                check(n, &kq, x, y);
            } else printf("Vui long tinh he so a, b truoc khi kiem tra ket qua!\n");
            break;
        case 4:
            if (t_check) {
            saveFile(n, kq, x, y);
            break;
            } else printf("Vui long tinh he so a, b truoc khi luu ket qua vao file!\n");
             break;
        case 5:
            printf("Quay lai menu chinh.\n");
            return;
        default:
            printf("Lua chon khong hop le. Vui long chon lai.\n");
        }
    }
    while (luachon2 != 5);
}

void menu() {
    Danhsachdac x, y;
    ketqua kq;
    char nameFile[256];
    int luachon, n = 0;

    do {
        printf("\n====================================== MENU ======================================\n");
        printf("1. Nhap du lieu tu file.\n");
        printf("2. Dung du lieu mau.\n");
        printf("3. Thoat.\n");
        printf("Lua chon: ");
        scanf("%d", &luachon);

        switch (luachon) {
            case 1:
                printf("Nhap du lieu tu file.\n");
                nhaptenfile(nameFile);
                n = readFile(nameFile, &x, &y);
                if (n<=0){
                    printf("Doc du lieu tu file that bai. Vui long kiem tra lai file va thu lai.\n");
                    break;
                }
                printf("Da doc duoc %d cap du lieu tu file %s.\n", n, nameFile);\
                menuin(luachon, x, y, kq, n);
                break;
            case 2:
                printf("Dung du lieu mau.\n");
                double mauX[] = {1.0, 2.0, 3.0, 4.0, 5.0};
                double mauY[] = {5.0, 320.0, 3645.0, 20480.0, 78125.0};
                n = 5;
                x.n = n;
                y.n = n;
                for (int i = 0; i < n; i++) {
                    x.data[i] = mauX[i];
                    y.data[i] = mauY[i];
                }
                printf("Da nap %d cap du lieu mau.\n", n);
                menuin(luachon, x, y, kq, n);
                break;
            case 3:
                printf("Tam biet!\n");
                break;
            default:
                printf("Lua chon khong hop le. Vui long chon lai.\n");
        }
    } while (luachon != 3);
}
