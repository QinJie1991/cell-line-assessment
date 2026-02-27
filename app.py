# -*- coding: utf-8 -*-
"""
慢病毒与细胞系构建智能评估系统 V1
- 2000bp包装极限版 + 双分支文献检索
- 整合HPA人类蛋白表达数据（自动下载）
"""

import streamlit as st
import pandas as pd
import requests
from Bio import Entrez
from datetime import datetime
import time
import re
from typing import Dict, List, Optional, Tuple
import json
from bs4 import BeautifulSoup
from fake_useragent import UserAgent
import urllib.parse
from collections import Counter
import hashlib
import openai
import os
import io
import zipfile

# ==================== 配置与初始化 ====================
st.set_page_config(
    page_title="慢病毒与细胞系构建智能评估系统 V1",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 安全密钥配置
try:
    NCBI_EMAIL = st.secrets["NCBI_EMAIL"]
    NCBI_API_KEY = st.secrets.get("NCBI_API_KEY", "")
    APP_PASSWORD = st.secrets.get("APP_PASSWORD", "")
    BAILIAN_API_KEY = st.secrets.get("BAILIAN_API_KEY") or st.secrets.get("DASHSCOPE_API_KEY", "")
    AI_MODEL = st.secrets.get("AI_MODEL", "qwen-plus")
except Exception:
    st.error("⚠️ 请先配置 Secrets（NCBI_EMAIL 等）")
    st.stop()

# 密码保护
if APP_PASSWORD:
    if 'authenticated' not in st.session_state:
        st.session_state.authenticated = False
    if not st.session_state.authenticated:
        pwd = st.text_input("🔒 请输入访问密码", type="password")
        if pwd == APP_PASSWORD:
            st.session_state.authenticated = True
            st.rerun()
        elif pwd:
            st.error("密码错误")
        st.stop()

# 初始化 session state
if 'analysis_results' not in st.session_state:
    st.session_state.analysis_results = None
if 'search_history' not in st.session_state:
    st.session_state.search_history = []

# ==================== Addgene 爬取模块 ====================

class AddgeneScraper:
    """Addgene 质粒爬取器 - 精简版（仅匹配基因名）"""
    
    def __init__(self):
        self.base_url = "https://www.addgene.org"
        self.ua = UserAgent()
        self.session = requests.Session()
        self.session.headers.update({
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5',
        })
    
    @st.cache_data(ttl=86400, show_spinner=False)
    def search_plasmids(_self, gene_symbol: str, max_results: int = 5) -> List[Dict]:
        """搜索 Addgene 质粒 - 仅返回基因名匹配的结果"""
        try:
            query = urllib.parse.quote(f"{gene_symbol}")
            search_url = f"{_self.base_url}/search/?q={query}&type=plasmid"
            headers = {'User-Agent': _self.ua.random}
            time.sleep(1)
            
            response = _self.session.get(search_url, headers=headers, timeout=15)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'lxml')
            plasmids = []
            
            result_items = (soup.find_all('article', class_='addgene-search-result') or 
                           soup.find_all('div', class_='search-result-item') or
                           soup.select('.plasmid-item'))
            
            if not result_items:
                result_items = soup.select('[data-testid="plasmid-card"]')
            
            for item in result_items[:max_results]:
                try:
                    plasmid_data = _self._parse_plasmid_card(item, gene_symbol)
                    if plasmid_data:
                        plasmids.append(plasmid_data)
                except Exception:
                    continue
            
            if not plasmids:
                plasmids = _self._search_by_gene_page(gene_symbol)
            
            return plasmids
            
        except Exception as e:
            st.error(f"Addgene 爬取错误: {str(e)}")
            return []
    
    def _parse_plasmid_card(self, card, gene_symbol: str) -> Optional[Dict]:
        """解析单个质粒卡片 - 仅提取基本信息和基因名匹配"""
        try:
            link_tag = card.find('a', href=re.compile(r'/\d{5,6}/')) or card.find('a', href=True)
            if not link_tag:
                return None
            
            href = link_tag.get('href', '')
            plasmid_id_match = re.search(r'/(\d{5,6})/', href)
            if not plasmid_id_match:
                return None
            
            plasmid_id = plasmid_id_match.group(1)
            full_url = f"{self.base_url}{href}" if href.startswith('/') else href
            
            name_tag = card.find('h3') or card.find('h2') or card.find('a', class_='title')
            name = name_tag.get_text(strip=True) if name_tag else "Unknown"
            
            # 严格匹配：基因名必须出现在质粒名称或描述中
            name_lower = name.lower()
            desc = card.get_text().lower()
            gene_lower = gene_symbol.lower()
            
            if gene_lower not in name_lower and gene_lower not in desc:
                return None
            
            return {
                'plasmid_id': plasmid_id,
                'name': name,
                'url': full_url,
                'insert_gene': gene_symbol,
            }
            
        except Exception:
            return None
    
    def _search_by_gene_page(self, gene_symbol: str) -> List[Dict]:
        """通过基因专属页面搜索"""
        try:
            gene_url = f"{self.base_url}/browse/gene/{gene_symbol}/"
            headers = {'User-Agent': self.ua.random}
            response = self.session.get(gene_url, headers=headers, timeout=10)
            
            if response.status_code != 200:
                return []
            
            soup = BeautifulSoup(response.text, 'lxml')
            plasmids = []
            seen_ids = set()
            
            links = soup.find_all('a', href=re.compile(r'/\d{5,6}/'))
            for link in links[:5]:
                href = link.get('href', '')
                match = re.search(r'/(\d{5,6})/', href)
                if match:
                    pid = match.group(1)
                    if pid not in seen_ids:
                        seen_ids.add(pid)
                        name = link.get_text(strip=True) or f"{gene_symbol} related"
                        
                        if gene_symbol.lower() not in name.lower():
                            continue
                            
                        plasmids.append({
                            'plasmid_id': pid,
                            'name': name,
                            'url': f"{self.base_url}/{pid}/",
                            'insert_gene': gene_symbol,
                            'source': 'Gene page'
                        })
            return plasmids
        except Exception:
            return []


# ==================== HPA 基因表达数据模块（自动下载版） ====================

class HPAGeneData:
    """
    基于本地 proteinatlas.tsv 的人类蛋白表达数据查询
    首次使用自动从 Human Protein Atlas 官网下载
    """
    
    def __init__(self, tsv_path: str = "data/proteinatlas.tsv"):
        self.tsv_path = tsv_path
        self.df = None
        self.available_columns = []
        
        # 期望的列（按优先级）
        self.desired_columns = {
            'Gene': 'Gene',
            'Ensembl': 'Ensembl',
            'Uniprot': 'Uniprot',
            'Subcellular main location': 'Subcellular main location',
            'Reliability': 'Reliability',
            'RNA tissue specific nTPM': 'RNA tissue specific nTPM'
        }
        
        # 如果文件不存在，尝试自动下载
        if not os.path.exists(self.tsv_path):
            self._auto_download()
        
        self._load_data()
    
    def _auto_download(self):
        """自动下载 HPA 数据文件（首次使用）"""
        try:
            st.warning("⬇️ 首次使用，正在下载 HPA 数据文件（约 30MB，请耐心等待）...")
            
            os.makedirs("data", exist_ok=True)
            url = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
            
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            response = requests.get(url, stream=True, timeout=300)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0
            chunk_size = 1024 * 1024
            
            zip_buffer = io.BytesIO()
            
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    zip_buffer.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        progress = min(downloaded / total_size, 1.0)
                        progress_bar.progress(progress)
                        status_text.text(f"下载进度: {downloaded/1024/1024:.1f} MB / {total_size/1024/1024:.1f} MB")
            
            status_text.text("正在解压...")
            
            zip_buffer.seek(0)
            with zipfile.ZipFile(zip_buffer, 'r') as zip_ref:
                zip_ref.extractall("data")
            
            progress_bar.empty()
            status_text.empty()
            
            if os.path.exists(self.tsv_path):
                st.success("✅ HPA 数据下载完成！")
                time.sleep(1)
            else:
                raise FileNotFoundError("解压后未找到 proteinatlas.tsv")
                
        except Exception as e:
            st.error(f"❌ 自动下载失败: {e}")
    
    def _load_data(self):
        """加载 TSV 文件（自动检测列名）"""
        try:
            if not os.path.exists(self.tsv_path):
                st.warning("⚠️ HPA 数据文件未找到")
                self.df = pd.DataFrame()
                return
            
            # 先读取第一行看看有哪些列
            sample_df = pd.read_csv(self.tsv_path, sep='\t', nrows=2)
            actual_columns = sample_df.columns.tolist()
            
            # 找出实际存在的列（不区分大小写，容忍空格差异）
            column_mapping = {}
            for desired_col in self.desired_columns.keys():
                # 精确匹配
                if desired_col in actual_columns:
                    column_mapping[desired_col] = desired_col
                else:
                    # 尝试不区分大小写匹配
                    for actual_col in actual_columns:
                        if actual_col.lower().strip() == desired_col.lower().strip():
                            column_mapping[desired_col] = actual_col
                            break
            
            if not column_mapping:
                st.error("HPA 文件中没有找到任何期望的列")
                self.df = pd.DataFrame()
                return
            
            # 只读取存在的列
            cols_to_read = list(column_mapping.values())
            self.df = pd.read_csv(
                self.tsv_path, 
                sep='\t', 
                low_memory=False,
                usecols=cols_to_read,
                dtype=str  # 全部作为字符串读取，避免类型问题
            )
            
            # 重命名为标准名称
            reverse_mapping = {v: k for k, v in column_mapping.items()}
            self.df = self.df.rename(columns=reverse_mapping)
            
            self.available_columns = list(column_mapping.keys())
            
            st.success(f"✅ HPA 数据已加载: {len(self.df):,} 条基因 (可用字段: {', '.join(self.available_columns)})")
            
        except Exception as e:
            st.error(f"加载 HPA 数据失败: {e}")
            self.df = pd.DataFrame()
    
    def get_gene_data(self, gene_symbol: str) -> Dict:
        """获取指定基因的表达数据"""
        if self.df is None or self.df.empty:
            return {}
        
        try:
            matches = self.df[self.df['Gene'].str.upper() == gene_symbol.upper()]
            
            if matches.empty:
                return {}
            
            row = matches.iloc[0]
            
            # 只返回实际存在的字段
            result = {}
            if 'Ensembl' in self.available_columns:
                result['ensembl_id'] = str(row.get('Ensembl', ''))
            if 'Uniprot' in self.available_columns:
                result['uniprot_id'] = str(row.get('Uniprot', ''))
            if 'Subcellular main location' in self.available_columns:
                result['subcellular_location'] = str(row.get('Subcellular main location', ''))
            if 'Reliability' in self.available_columns:
                result['reliability'] = str(row.get('Reliability', ''))
            if 'RNA tissue specific nTPM' in self.available_columns:
                result['rna_tissue_specificity'] = str(row.get('RNA tissue specific nTPM', ''))
            
            # 构建链接（如果有 Ensembl）
            if 'ensembl_id' in result and result['ensembl_id']:
                result['hpa_link'] = f"https://www.proteinatlas.org/{result['ensembl_id']}"
            
            return result
            
        except Exception as e:
            st.error(f"查询 HPA 基因数据失败: {e}")
            return {}
    
    def check_data_available(self) -> bool:
        """检查数据是否可用"""
        return self.df is not None and not self.df.empty


# ==================== NCBI 数据获取模块 ====================

class BioDataFetcher:
    def __init__(self, email: str, api_key: str = ""):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.uniprot_base = "https://rest.uniprot.org/uniprotkb/search.json"
        self.headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'}
        self.addgene_scraper = AddgeneScraper()
        # 关键：使用 HPAGeneData（不是 HumanAtlasScraper）
        self.hpa_data = HPAGeneData()
    
    def get_ncbi_gene_info(self, gene_symbol: str, species: str) -> Dict:
        try:
            term = f"{gene_symbol}[Gene Name] AND {species}[Organism]"
            handle = Entrez.esearch(db="gene", term=term, retmax=1)
            record = Entrez.read(handle)
            handle.close()
            
            if not record["IdList"]:
                return {"status": "not_found", "error": f"未找到 {gene_symbol} ({species})"}
            
            gene_id = record["IdList"][0]
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
            gene_data = Entrez.read(handle)
            handle.close()
            
            gene_entry = gene_data[0]
            summary = gene_entry.get("Entrezgene_summary", "")
            
            lethal_keywords = ["essential", "lethal", "required for cell viability", 
                             "knockout mice die", "embryonic lethal"]
            phenotype = "非必需"
            if any(kw in summary.lower() for kw in lethal_keywords):
                phenotype = "必需（潜在致死风险）"
            
            return {
                "gene_id": gene_id,
                "symbol": gene_symbol,
                "species": species,
                "description": summary[:500] if summary else "无描述",
                "phenotype": phenotype,
                "chromosome": gene_entry.get("Entrezgene_location", [{}])[0].get("Gene-location", {}).get("Gene-location_chromosome", "N/A"),
                "status": "success"
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}
    
    def get_uniprot_info(self, gene_symbol: str, species: str) -> Dict:
        try:
            species_map = {"Homo sapiens": "human", "Mus musculus": "mouse", "Rattus norvegicus": "rat"}
            org = species_map.get(species, species.lower())
            
            query = f"gene:{gene_symbol}+organism:{org}"
            params = {
                "query": query,
                "fields": "accession,gene_names,length,cc_subcellular_location,sequence",
                "format": "json",
                "size": 1
            }
            
            response = requests.get(self.uniprot_base, params=params, headers=self.headers, timeout=10)
            data = response.json()
            
            if not data.get("results"):
                return {"status": "not_found", "error": f"UniProt 未找到 {gene_symbol}"}
            
            protein = data["results"][0]
            accession = protein.get("primaryAccession", "")
            seq_length = protein.get("sequence", {}).get("length", 0)
            
            loc_text = ""
            comments = protein.get("comments", [])
            for comment in comments:
                if comment.get("commentType") == "SUBCELLULAR LOCATION":
                    locations = comment.get("subcellularLocations", [])
                    locs = [loc.get("location", {}).get("value", "") for loc in locations]
                    loc_text = "; ".join([l for l in locs if l])
            
            cds_length = seq_length * 3 if seq_length else 0
            
            return {
                "uniprot_id": accession,
                "protein_length": seq_length,
                "cds_length_bp": cds_length,
                "subcellular_location": loc_text or "未标注",
                "status": "success"
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}


# ==================== 通义千问 AI 分析模块 ====================

class AIAnalyzer:
    def __init__(self):
        self.client = None
        self.model = None
        self.cache = {}
        
        try:
            api_key = BAILIAN_API_KEY
            
            if not api_key:
                raise ValueError("未配置 BAILIAN_API_KEY 或 DASHSCOPE_API_KEY")
            
            self.model = AI_MODEL
            
            self.client = openai.OpenAI(
                api_key=api_key,
                base_url="https://dashscope.aliyuncs.com/compatible-mode/v1"
            )
            
            st.success(f"✅ 通义千问已连接（模型：{self.model}）")
            
        except Exception as e:
            st.error(f"通义千问初始化失败：{e}")
            raise e
    
    def _get_cache_key(self, text: str, task: str) -> str:
        return hashlib.md5(f"{task}:{text}".encode()).hexdigest()[:16]
    
    def _call_qwen(self, prompt: str, json_mode: bool = True, temperature: float = 0.3) -> Dict:
        try:
            messages = [
                {"role": "system", "content": "你是资深细胞生物学和分子克隆专家，擅长分析基因功能、设计细胞系构建实验方案。请用中文回答，专业且具体。"},
                {"role": "user", "content": prompt}
            ]
            
            kwargs = {
                "model": self.model,
                "messages": messages,
                "temperature": temperature,
                "max_tokens": 2000
            }
            
            if json_mode:
                kwargs["response_format"] = {"type": "json_object"}
            
            with st.spinner("🤖 通义千问分析中..."):
                response = self.client.chat.completions.create(**kwargs)
                content = response.choices[0].message.content
                
                if json_mode:
                    return json.loads(content)
                return {"result": content}
                
        except Exception as e:
            st.error(f"通义千问API调用失败：{str(e)}")
            return {"error": str(e)}
    
    def assess_gene_essentiality(self, gene_data: Dict, uniprot_data: Dict) -> Dict:
        gene_name = gene_data.get("symbol", "")
        description = gene_data.get("description", "")
        location = uniprot_data.get("subcellular_location", "")
        
        cache_key = self._get_cache_key(f"{gene_name}_{description[:100]}", "essentiality")
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        prompt = f"""请基于以下基因信息，分析该基因是否为细胞必需基因（essential gene），并评估构建KO细胞系的可行性。

基因名称：{gene_name}
基因功能描述：{description}
亚细胞定位：{location}

请按以下JSON格式回答：
{{
    "is_essential": "是/否/可能",
    "confidence": "高/中/低",
    "rationale": "详细解释判断理由，基于该基因的生物学功能",
    "expected_phenotype": "如果敲除，预期会出现什么细胞表型（如细胞死亡、增殖减慢、周期阻滞等）",
    "construction_difficulty": "容易/中等/困难",
    "recommendation": "具体的实验建议（如是否建议使用诱导型系统、先尝试KD还是直接KO）",
    "key_considerations": "实验设计中的关键注意事项（2-3点）",
    "alternative_approaches": "如果KO困难，建议的替代方案（如使用siRNA、诱导型敲除等）"
}}

注意：
1. 如果基因是管家基因（如GAPDH、ACTB）或参与DNA复制/修复核心功能，通常判定为必需
2. 如果是癌基因或信号通路分子，通常为非必需但可能有表型
3. 请给出具体生物学机制解释，不要泛泛而谈"""
        
        result = self._call_qwen(prompt, json_mode=True)
        self.cache[cache_key] = result
        return result
    
    def analyze_literature_deep(self, articles: List[Dict], gene: str, construct_type: str) -> Dict:
        if not articles:
            return {"summary": "无相关文献", "protocols": []}
        
        articles_text = "\n\n".join([
            f"文献{i+1}:\n标题：{a['title']}\n方法关键词：{a['methods']}\n摘要片段：{a.get('abstract_snippet', '')}"
            for i, a in enumerate(articles[:5])
        ])
        
        cache_key = self._get_cache_key(articles_text[:500], f"lit_{construct_type}")
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        prompt = f"""作为方法学专家，分析以下关于"{gene}"基因{construct_type}（过表达/敲低/敲除）细胞系构建的文献。

{articles_text}

请提取以下信息并以JSON格式返回：
{{
    "summary": "这些研究的总体方法学趋势（中文，2-3句话概括）",
    "protocols": [
        {{
            "cell_line": "使用的细胞系",
            "vector_system": "具体载体名称（如pLKO.1, lentiCRISPRv2, pLVX等）",
            "delivery_method": "递送方法（如慢病毒转导、脂质体转染Lipofectamine 2000/3000、电转等）",
            "selection_method": "筛选方法和浓度（如Puromycin 2μg/mL, G418 500μg/mL等）",
            "validation_method": "验证方法（如Western Blot、qPCR、免疫荧光、流式等）",
            "efficiency": "效率信息（如有提及）",
            "duration": "实验周期（如有提及）"
        }}
    ],
    "common_methods": "最常见的方法组合（中文总结）",
    "critical_factors": "影响实验成功的关键因素（2-3点）",
    "troubleshooting": "常见问题及解决方案（如污染、效率低、细胞死亡等）",
    "recommendation": "针对该基因{construct_type}的具体操作建议"
}}

注意：
1. 如果某字段在文献中未提及，填写"未提及"
2. 优先提取具体的试剂名称和浓度（如Lipofectamine 2000, Polybrene 8μg/mL等）
3. 如果提到多种细胞系，列出最常用的2-3种"""
        
        result = self._call_qwen(prompt, json_mode=True)
        self.cache[cache_key] = result
        return result


# ==================== 分析主类 ====================

class ConstructAnalyzer:
    def __init__(self):
        self.fetcher = BioDataFetcher(NCBI_EMAIL, NCBI_API_KEY)
    
    def analyze_gene(self, gene_symbol: str, species: str, 
                    cell_line: Optional[str] = None, 
                    cell_species: Optional[str] = None) -> Dict:
        
        with st.spinner(f"🔍 正在深度分析 {gene_symbol}..."):
            
            st.text("检索 NCBI Gene...")
            ncbi_info = self.fetcher.get_ncbi_gene_info(gene_symbol, species)
            time.sleep(0.5)
            
            st.text("检索 UniProt...")
            uniprot_info = self.fetcher.get_uniprot_info(gene_symbol, species)
            time.sleep(0.5)
            
            st.text("检索 Addgene...")
            addgene_plasmids = self.fetcher.addgene_scraper.search_plasmids(gene_symbol)
            time.sleep(0.5)
            
            # 仅人类基因查询HPA数据
            hpa_gene_data = {}
            if species == "Homo sapiens":
                st.text("检索 HPA 蛋白表达数据...")
                hpa_gene_data = self.fetcher.hpa_data.get_gene_data(gene_symbol)
            time.sleep(0.5)
            
            lentiviral = self._assess_lentiviral(ncbi_info, uniprot_info, addgene_plasmids)
            literature = self._search_all_constructs(gene_symbol, cell_line)
            
            ai_analysis = {}
            try:
                ai_analyzer = AIAnalyzer()
                
                with st.spinner("🤖 通义千问正在深度分析..."):
                    ai_analysis["gene_assessment"] = ai_analyzer.assess_gene_essentiality(
                        ncbi_info, uniprot_info
                    )
                    
                    # AI分析特定细胞的文献（如果有）
                    if cell_line and literature.get("specific_cell", {}).get("found"):
                        for ctype in ["overexpression", "knockdown", "knockout"]:
                            if literature["specific_cell"][ctype]["articles"]:
                                ai_result = ai_analyzer.analyze_literature_deep(
                                    literature["specific_cell"][ctype]["articles"], gene_symbol, ctype
                                )
                                literature["specific_cell"][ctype]["ai_analysis"] = ai_result
                    
                    # AI分析通用文献
                    for ctype in ["overexpression", "knockdown", "knockout"]:
                        if literature[ctype]["articles"]:
                            ai_result = ai_analyzer.analyze_literature_deep(
                                literature[ctype]["articles"], gene_symbol, ctype
                            )
                            literature[ctype]["ai_analysis"] = ai_result
                    
                    ai_analysis["enabled"] = True
                    
            except Exception as e:
                st.warning(f"AI分析未启用：{e}")
                ai_analysis["enabled"] = False
            
            result = {
                "input_info": {
                    "gene_symbol": gene_symbol,
                    "species": species,
                    "cell_line": cell_line or "未指定",
                    "cell_species": cell_species or "未指定",
                    "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M"),
                    "ai_enabled": ai_analysis.get("enabled", False)
                },
                "gene_function": ncbi_info,
                "protein_data": uniprot_info,
                "hpa_gene_data": hpa_gene_data,
                "addgene_plasmids": addgene_plasmids,
                "lentiviral_assessment": lentiviral,
                "cell_line_constructs": literature,
                "ai_analysis": ai_analysis,
                "database_record": self._format_database_record(
                    gene_symbol, species, cell_line, ncbi_info, uniprot_info, 
                    lentiviral, literature, addgene_plasmids, hpa_gene_data
                )
            }
            
            return result
    
    def _assess_lentiviral(self, ncbi_info: Dict, uniprot_info: Dict, plasmids: List) -> Dict:
        """评估慢病毒适用性 - 2000bp极限简化版"""
        warnings = []
        recommendations = []
        cds_len = uniprot_info.get("cds_length_bp", 0)
        
        # 致死性判断
        if ncbi_info.get("phenotype") == "必需（潜在致死风险）":
            warnings.append("🚨 必需基因，敲除可能导致细胞致死，建议使用诱导型系统")
        
        # 长度判断
        if cds_len > 2000:
            warnings.append(f"⚠️ 长度超限（{cds_len}bp > 2000bp）")
            warnings.append("💡 建议：去除荧光标签或使用split-vector系统")
            suitable = False
            rating = "❌ 不推荐（长度超限）"
        else:
            suitable = True
            rating = "✅ 可接受"
        
        # Addgene资源提示
        if plasmids:
            recommendations.append(f"Addgene 提供 {len(plasmids)} 个质粒")
        
        return {
            "suitable": suitable,
            "cds_length": cds_len,
            "packaging_limit": 2000,
            "warnings": warnings,
            "recommendations": recommendations,
            "overall_assessment": rating
        }
    
    def _search_all_constructs(self, gene_symbol: str, cell_line: Optional[str]) -> Dict:
        """双分支文献检索"""
        results = {}
        
        # 分支1：特定目的细胞系构建
        if cell_line:
            query_specific = f'{gene_symbol}[Title/Abstract] AND {cell_line}[Title/Abstract] AND (cell line OR cell-line)'
            
            try:
                handle = Entrez.esearch(db="pubmed", term=query_specific, retmax=100, sort="relevance")
                record = Entrez.read(handle)
                handle.close()
                
                pmids = record["IdList"]
                
                if pmids:
                    fetch_ids = pmids[:20]
                    handle = Entrez.efetch(db="pubmed", id=fetch_ids, rettype="abstract", retmode="xml")
                    articles = Entrez.read(handle)
                    handle.close()
                    
                    specific_oe, specific_kd, specific_ko = [], [], []
                    
                    for article in articles.get("PubmedArticle", []):
                        try:
                            medline = article["MedlineCitation"]
                            article_data = medline["Article"]
                            title = article_data.get("ArticleTitle", "").lower()
                            abstract = str(article_data.get("Abstract", {}).get("AbstractText", "")).lower()
                            full_text = title + " " + abstract
                            
                            pmid = str(medline.get("PMID", "N/A"))
                            title_display = article_data.get("ArticleTitle", "N/A")
                            
                            methods = [kw for kw in ["lentiviral", "transfection", "electroporation"] if kw in full_text]
                            
                            article_info = {
                                "pmid": pmid,
                                "title": title_display,
                                "methods": ", ".join(methods) if methods else "未明确提及",
                                "cell_line": cell_line
                            }
                            
                            if any(kw in full_text for kw in ["overexpression", "over-expression", "ectopic"]):
                                specific_oe.append(article_info)
                            elif any(kw in full_text for kw in ["knockdown", "sirna", "shrna"]):
                                specific_kd.append(article_info)
                            elif any(kw in full_text for kw in ["knockout", "crispr", "knock-out"]):
                                specific_ko.append(article_info)
                                    
                        except Exception:
                            continue
                    
                    results["specific_cell"] = {
                        "found": True,
                        "total_count": len(pmids),
                        "overexpression": {"articles": specific_oe[:5], "count": len(specific_oe)},
                        "knockdown": {"articles": specific_kd[:5], "count": len(specific_kd)},
                        "knockout": {"articles": specific_ko[:5], "count": len(specific_ko)},
                        "message": f"在 {cell_line} 中找到 {len(pmids)} 篇文献（OE:{len(specific_oe)} KD:{len(specific_kd)} KO:{len(specific_ko)}）"
                    }
                else:
                    results["specific_cell"] = {"found": False, "message": f"未在 {cell_line} 中找到相关研究"}
            except Exception as e:
                results["specific_cell"] = {"found": False, "message": f"检索失败: {e}"}
        else:
            results["specific_cell"] = {"found": False, "message": "未输入细胞系名称"}
        
        # 分支2：通用表达调控文献
        for construct_type in ["overexpression", "knockdown", "knockout"]:
            type_map = {
                "overexpression": "overexpression OR over-expression OR ectopic",
                "knockdown": "knockdown OR siRNA OR shRNA",
                "knockout": "knockout OR CRISPR OR knock-out"
            }
            
            query = f'{gene_symbol}[Title/Abstract] AND ({type_map[construct_type]}) AND (cell line OR cell-line)'
            
            try:
                handle = Entrez.esearch(db="pubmed", term=query, retmax=100, sort="relevance")
                record = Entrez.read(handle)
                handle.close()
                
                pmids = record["IdList"]
                articles_list = []
                
                if pmids:
                    fetch_ids = pmids[:10]
                    handle = Entrez.efetch(db="pubmed", id=fetch_ids, rettype="abstract", retmode="xml")
                    articles = Entrez.read(handle)
                    handle.close()
                    
                    for article in articles.get("PubmedArticle", []):
                        try:
                            medline = article["MedlineCitation"]
                            article_data = medline["Article"]
                            title = article_data.get("ArticleTitle", "N/A")
                            pmid = str(medline.get("PMID", "N/A"))
                            abstract = str(article_data.get("Abstract", {}).get("AbstractText", ""))
                            
                            methods = [kw for kw in ["lentiviral", "transfection", "electroporation"] 
                                      if kw in (title + abstract).lower()]
                            
                            articles_list.append({
                                "pmid": pmid,
                                "title": title,
                                "methods": ", ".join(methods) if methods else "未明确提及",
                                "abstract_snippet": abstract[:300] + "..." if len(abstract) > 300 else abstract
                            })
                        except Exception:
                            continue
                
                results[construct_type] = {"count": len(pmids), "articles": articles_list}
                
            except Exception as e:
                results[construct_type] = {"count": 0, "articles": [], "error": str(e)}
        
        return results
    
    def _format_database_record(self, gene_symbol: str, species: str, cell_line: Optional[str],
                               ncbi_info: Dict, uniprot_info: Dict, lentiviral: Dict,
                               literature: Dict, plasmids: List, hpa_gene_data: Dict) -> Dict:
        """格式化数据库记录"""
        return {
            "gene_symbol": gene_symbol,
            "species": species,
            "cell_line": cell_line,
            "gene_id": ncbi_info.get("gene_id"),
            "uniprot_id": uniprot_info.get("uniprot_id"),
            "ensembl_id": hpa_gene_data.get("ensembl_id"),
            "cds_length": uniprot_info.get("cds_length_bp"),
            "is_essential": ncbi_info.get("phenotype") == "必需（潜在致死风险）",
            "lentiviral_suitable": lentiviral.get("suitable"),
            "literature_count": {
                "oe": literature.get("overexpression", {}).get("count", 0),
                "kd": literature.get("knockdown", {}).get("count", 0),
                "ko": literature.get("knockout", {}).get("count", 0)
            },
            "plasmid_count": len(plasmids),
            "hpa_data_available": bool(hpa_gene_data),
            "timestamp": datetime.now().isoformat()
        }


# ==================== Streamlit UI ====================

def main():
    st.title("🔬 慢病毒与细胞系构建智能评估系统 V1")
    st.markdown("**2000bp包装极限版 + 双分支文献检索 + AI 智能分析**")
    
    with st.sidebar:
        st.header("⚙️ 分析参数设置")
        
        gene_symbol = st.text_input("基因符号 (Gene Symbol)", "TP53").strip().upper()
        
        species = st.selectbox(
            "基因来源物种",
            ["Homo sapiens", "Mus musculus", "Rattus norvegicus"],
            index=0
        )
        
        cell_line = st.text_input("目的细胞系 (可选)", "HEK293").strip()
        
        cell_species = st.selectbox(
            "细胞系物种",
            ["未指定", "Human", "Mouse", "Rat", "Other"],
            index=0
        )
        
        analyze_btn = st.button("🚀 开始深度分析", type="primary", use_container_width=True)
        
        st.divider()
        
        # HPA 数据状态提示
        hpa_checker = HPAGeneData()  # 临时实例检查状态
        if hpa_checker.check_data_available():
            st.success("✅ HPA 人类蛋白数据已加载")
        else:
            st.warning("⚠️ HPA 数据未配置")
            st.caption("首次使用将自动下载（约30MB）")
        
        st.info("""
        **系统功能：**
        - ✅ NCBI/UniProt 基因信息检索
        - ✅ Addgene 质粒资源查询
        - ✅ HPA 人类蛋白表达数据（自动下载）
        - ✅ 2000bp 慢病毒包装极限评估
        - ✅ PubMed 双分支文献检索
        - ✅ 通义千问 AI 智能分析
        """)
    
    if analyze_btn and gene_symbol:
        analyzer = ConstructAnalyzer()
        
        try:
            result = analyzer.analyze_gene(
                gene_symbol=gene_symbol,
                species=species,
                cell_line=cell_line if cell_line else None,
                cell_species=cell_species if cell_species != "未指定" else None
            )
            
            st.session_state.analysis_results = result
            st.session_state.search_history.append({
                "gene": gene_symbol,
                "cell": cell_line,
                "time": datetime.now().strftime("%H:%M")
            })
            
        except Exception as e:
            st.error(f"分析过程出错: {e}")
    
    # 显示结果
    if st.session_state.analysis_results:
        display_results(st.session_state.analysis_results)

def display_results(result: Dict):
    """展示分析结果"""
    info = result["input_info"]
    
    st.header(f"📊 {info['gene_symbol']} 分析报告")
    cols = st.columns(4)
    cols[0].metric("基因", info['gene_symbol'])
    cols[1].metric("物种", info['species'].split()[0])
    cols[2].metric("细胞系", info['cell_line'])
    cols[3].metric("AI分析", "已启用" if info['ai_enabled'] else "未启用")
    
    tabs = st.tabs(["🧬 基因功能", "🦠 慢病毒评估", "📚 文献检索", "🧪 实验资源"])
    
    # Tab 1: 基因功能（包含HPA数据）
    with tabs[0]:
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.subheader("NCBI Gene")
            gene_data = result["gene_function"]
            if gene_data.get("status") == "success":
                st.write(f"**基因ID:** {gene_data.get('gene_id')}")
                st.write(f"**染色体:** {gene_data.get('chromosome')}")
                st.write(f"**必需性:** {gene_data.get('phenotype')}")
                with st.expander("功能描述"):
                    st.write(gene_data.get("description", "无"))
            else:
                st.error(gene_data.get("error", "无法获取基因信息"))
        
        with col2:
            st.subheader("UniProt")
            prot_data = result["protein_data"]
            if prot_data.get("status") == "success":
                st.write(f"**UniProt ID:** {prot_data.get('uniprot_id')}")
                st.write(f"**蛋白长度:** {prot_data.get('protein_length')} aa")
                st.write(f"**CDS长度:** {prot_data.get('cds_length_bp')} bp")
                with st.expander("亚细胞定位"):
                    st.write(prot_data.get("subcellular_location", "未标注"))
            else:
                st.error(prot_data.get("error", "无法获取蛋白信息"))
        
        with col3:
            st.subheader("HPA 蛋白表达")
            hpa_data = result.get("hpa_gene_data", {})
            if hpa_data:
                st.write(f"**Ensembl:** {hpa_data.get('ensembl_id', 'N/A')}")
                st.write(f"**可靠性:** {hpa_data.get('reliability', 'N/A')}")
                rna_expr = hpa_data.get('rna_tissue_specificity', 'N/A')
                if len(str(rna_expr)) > 50:
                    rna_expr = str(rna_expr)[:47] + "..."
                st.write(f"**RNA表达:** {rna_expr}")
                with st.expander("亚细胞定位"):
                    st.write(hpa_data.get('subcellular_location', '未标注'))
                st.caption(f"[HPA详情页]({hpa_data.get('hpa_link', '#')})")
            else:
                if info['species'] == "Homo sapiens":
                    st.info("未找到HPA数据")
                else:
                    st.info("HPA仅支持人类基因")
    
    # Tab 2: 慢病毒评估
    with tabs[1]:
        lentiviral = result["lentiviral_assessment"]
        
        col1, col2, col3 = st.columns(3)
        col1.metric("CDS长度", f"{lentiviral['cds_length']} bp")
        col2.metric("包装极限", "2000 bp")
        col3.metric("评估结果", lentiviral['overall_assessment'])
        
        if lentiviral['warnings']:
            st.error("**警告:**")
            for warning in lentiviral['warnings']:
                st.write(warning)
        
        if lentiviral['recommendations']:
            st.info("**建议:**")
            for rec in lentiviral['recommendations']:
                st.write(f"- {rec}")
        
        st.divider()
        st.markdown("""
        **2000bp 极限说明:**
        - 第三代慢病毒系统包装容量约 8kb（含载体骨架）
        - 实际插入片段建议 ≤2000bp 以确保滴度
        - 超过 2000bp 建议：使用 split-vector 系统或选择其他递送方式
        """)
    
    # Tab 3: 文献检索
    with tabs[2]:
        literature = result["cell_line_constructs"]
        
        if literature.get("specific_cell", {}).get("found"):
            st.success(literature["specific_cell"]["message"])
            
            subtabs = st.tabs(["过表达", "敲低", "敲除"])
            types_map = ["overexpression", "knockdown", "knockout"]
            
            for tab, ctype in zip(subtabs, types_map):
                with tab:
                    data = literature["specific_cell"][ctype]
                    st.write(f"找到 {data['count']} 篇相关文献")
                    
                    if data["articles"]:
                        for article in data["articles"]:
                            with st.expander(f"{article['title'][:80]}..."):
                                st.write(f"**PMID:** {article['pmid']}")
                                st.write(f"**方法:** {article['methods']}")
                    
                    if "ai_analysis" in data:
                        st.markdown("**🤖 AI 方法学分析**")
                        ai_data = data["ai_analysis"]
                        st.write(ai_data.get("summary", ""))
                        if ai_data.get("common_methods"):
                            st.info(f"常用方法: {ai_data['common_methods']}")
        else:
            st.warning(literature.get("specific_cell", {}).get("message", "未检索特定细胞系文献"))
        
        st.divider()
        st.subheader("通用基因表达调控文献")
        
        cols = st.columns(3)
        for idx, ctype in enumerate(["overexpression", "knockdown", "knockout"]):
            with cols[idx]:
                count = literature[ctype]["count"]
                label = ctype.replace("overexpression", "OE").replace("knockdown", "KD").replace("knockout", "KO")
                st.metric(label, f"{count} 篇")
    
    # Tab 4: 实验资源
    with tabs[3]:
        st.subheader("Addgene 质粒资源")
        plasmids = result["addgene_plasmids"]
        
        if plasmids:
            for p in plasmids:
                st.write(f"**[{p['plasmid_id']}]** {p['name']}")
                st.caption(f"[查看详情]({p['url']})")
                st.divider()
        else:
            st.info("未找到相关质粒")

if __name__ == "__main__":
    main()

