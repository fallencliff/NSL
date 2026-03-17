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

    %% =2: 中心控制层 (Master Node) =
    subgra
