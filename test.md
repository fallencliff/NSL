```mermaid
flowchart TB
    %% ================= 样式定义 =================
    classDef userLayer fill:#f9f9f9,stroke:#333,stroke-width:2px;
    classDef masterZone fill:#e3f2fd,stroke:#1e88e5,stroke-width:2px;
    classDef mqLayer fill:#fff3e0,stroke:#fb8c00,stroke-width:2px;
    classDef probeZone fill:#e8f5e9,stroke:#43a047,stroke-width:2px;
    classDef targetZone fill:#fce4ec,stroke:#e53935,stroke-width:2px;
    classDef firewall fill:#ffebee,stroke:#c62828,stroke-width:4px,stroke-dasharray: 5 5;
    classDef external fill:#ede7f6,stroke:#5e35b1,stroke-width:2px,stroke-dasharray: 5 5;

    %% ================= L1: 用户与外部依赖层 =================
    subgraph Users_External ["👨‍💻 访问层 & 外部系统"]
        U1(运维主管/总监)
        U2(一线运维/专家)
        SSO[[🏢 4A SSO 认证系统]]:::external
        LLM[[🧠 内部私有大模型 API]]:::external
        MSG[[📧 邮件/内部社交工具推送]]:::external
        U1 & U2 -. 登录鉴权 .-> SSO
    end

    %% ================= L2: 中心控制层 (Master Node) =================
    subgraph Master_Control ["🏢 核心管理区 (Master Node - Python 体系)"]
        direction TB
        Gateway{API 网关 / Nginx<br/>TLS, 限流, 鉴权}
        
        subgraph Master_Services ["核心微服务群 (FastAPI / Celery)"]
            UI[🖥️ 全景驾驶舱前端<br/>Vue3/大屏/像素矩阵]
            CMDB_Router[🗃️ 资产与探针路由引擎<br/>Zone/Probe映射]
            Scheduler[⏱️ 动态调度引擎<br/>Celery/重保日历Cron]
            Vault[🔒 凭证金库 Vault<br/>AES-256 加密/RSA 二加密]
            AI_Pipeline[🛡️ AI 数据管道<br/>正则脱敏/Prompt组装]
        end

        Gateway --> UI
        Gateway --> CMDB_Router
        Gateway --> Scheduler
        Gateway --> Vault
        Gateway --> AI_Pipeline
    end

    %% ================= L3: 数据与通信总线层 =================
    subgraph Data_Bus ["🚄 通信与数据底座"]
        direction LR
        MQ((🐰 RabbitMQ<br/>指令分发/异步队列)):::mqLayer
        gRPC_S((📡 gRPC Server<br/>探针心跳/文件回传)):::mqLayer
        DB_PostgreSQL[(PostgreSQL<br/>台账/元数据)]
        DB_Elasticsearch[(Elasticsearch<br/>历史快照/JSON日志)]
        DB_Redis[(Redis<br/>状态缓存/分布式锁)]
    end

    %% ================= L4: 严格的网络隔离边界 =================
    FW==========[🚨 核心防火墙策略：仅允许各区探针主动(反向)连接 MQ/gRPC 端口]==========:::firewall

    %% ================= L5: 边缘探针层 (Probe Nodes) =================
    subgraph Edge_Probes ["🌐 边缘探针层 (GoLang 轻量探针) - 部署于各区跳板机"]
        direction LR
        Probe_A[探针 A<br/>(省公司业务网)]
        Probe_B[探针 B<br/>(省公司BMC带外网)]
        Probe_C[探针 C<br/>(集团一级IT云)]
    end

    %% ================= L6: 目标执行层 (Agentless Targets) =================
    subgraph Targets_Exec ["🎯 无代理目标层"]
        direction LR
        Target_OS[💻 业务网服务器<br/>(SSH/JDBC)]
        Target_BMC[🎛️ 物理机硬件BMC<br/>(IPMI/Redfish)]
        Target_Cloud[☁️ 集团云虚机/裸金属<br/>(SSH)]
    end

    %% ================= 连线关系 (Flows) =================
    U1 & U2 == HTTPS ==> Gateway

    AI_Pipeline -. 请求脱敏数据 .-> LLM
    Scheduler -. 发送告警报告 .-> MSG

    Master_Services === MQ & gRPC_S & DB_PostgreSQL & DB_Elasticsearch & DB_Redis

    %% 跨防火墙核心链路 (探针主动反向连接)
    Probe_A & Probe_B & Probe_C == "TCP/TLS 主动连接\n(心跳/任务拉取/结果回传)" === FW
    FW === MQ & gRPC_S

    %% 探针直连目标 (内网就近管理)
    Probe_A -->|SSH/JDBC 并发直连| Target_OS
    Probe_B -->|IPMI/Redfish 并发直连| Target_BMC
    Probe_C -->|SSH/JDBC 并发直连| Target_Cloud

    %% ================= 样式应用 =================
    class U1,U2 userLayer;
    class Master_Control masterZone;
    class Users_External external;
    class Data_Bus mqLayer;
    class Edge_Probes probeZone;
    class Targets_Exec targetZone;
