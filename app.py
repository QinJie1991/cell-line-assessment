# -*- coding: utf-8 -*-
"""
慢病毒与细胞系构建智能评估系统 V1
- 2000bp包装极限版 + 双分支文献检索
- 整合HPA人类蛋白表达数据（自动下载）
- 增强慢病毒风险评估（功能+文献+序列）
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
    """Addgene 质粒爬取器 - 精简版"""
    
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
        """搜索 Addgene 质粒"""
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
        """解析单个质粒卡片"""
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


# ==================== HPA 基因表达数据模块 ====================

class HPAGeneData:
    """基于本地 proteinatlas.tsv 的人类蛋白表达数据查询"""
    
    def __init__(self, tsv_path: str = "data/proteinatlas.tsv"):
        self.tsv_path = tsv_path
        self.df = None
        self.available_columns = []
        self._load_data()
    
    def _load_data(self):
        """加载 TSV 文件（兼容多版本）"""
        try:
            if not os.path.exists(self.tsv_path):
                self._auto_download()
                if not os.path.exists(self.tsv_path):
                    return
            
            self.df = pd.read_csv(self.tsv_path, sep='\t', low_memory=False)
            self.available_columns = self.df.columns.tolist()
            self.df.columns = [col.strip() for col in self.df.columns]
            
            if 'Gene' not in self.df.columns:
                gene_col = None
                for col in self.df.columns:
                    if 'gene' in col.lower():
                        gene_col = col
                        break
                if gene_col:
                    self.df = self.df.rename(columns={gene_col: 'Gene'})
                else:
                    st.error("HPA 文件中找不到基因名列")
                    self.df = pd.DataFrame()
                    return
            
            st.success(f"✅ HPA 数据已加载: {len(self.df):,} 条基因")
            
        except Exception as e:
            st.error(f"加载 HPA 数据失败: {e}")
            self.df = pd.DataFrame()
    
    def _auto_download(self):
        """自动下载 HPA 数据文件"""
        try:
            st.info("⬇️ 正在下载 HPA 数据...")
            os.makedirs("data", exist_ok=True)
            
            url = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
            response = requests.get(url, timeout=300)
            
            with open("data/hpa_temp.zip", "wb") as f:
                f.write(response.content)
            
            with zipfile.ZipFile("data/hpa_temp.zip", 'r') as zip_ref:
                zip_ref.extractall("data")
            
            os.remove("data/hpa_temp.zip")
            st.success("✅ HPA 数据下载完成")
            
        except Exception as e:
            st.error(f"下载失败: {e}")
    
    def get_gene_data(self, gene_symbol: str) -> Dict:
        """获取基因数据"""
        if self.df is None or self.df.empty:
            return {}
        
        try:
            mask = self.df['Gene'].str.upper() == gene_symbol.upper()
            if not mask.any():
                return {}
            
            row = self.df[mask].iloc[0]
            result = {}
            
            if 'Ensembl' in self.df.columns:
                result['ensembl_id'] = str(row['Ensembl'])
            if 'Uniprot' in self.df.columns:
                result['uniprot_id'] = str(row['Uniprot'])
            if 'Subcellular main location' in self.df.columns:
                result['subcellular_location'] = str(row['Subcellular main location'])
            elif 'Subcellular location' in self.df.columns:
                result['subcellular_location'] = str(row['Subcellular location'])
            
            rna_col = None
            for col in self.df.columns:
                if 'tissue' in col.lower() and 'rna' in col.lower():
                    rna_col = col
                    break
            if rna_col:
                result['rna_tissue_specificity'] = str(row[rna_col])
            
            if 'Reliability' in self.df.columns:
                result['reliability'] = str(row['Reliability'])
            else:
                result['reliability'] = 'N/A'
            
            if 'ensembl_id' in result:
                result['hpa_link'] = f"https://www.proteinatlas.org/{result['ensembl_id']}"
            
            return result
            
        except Exception as e:
            st.error(f"查询 HPA 数据错误: {e}")
            return {}
    
    def check_data_available(self) -> bool:
        return self.df is not None and not self.df.empty


# ==================== 慢病毒风险评估类 ====================

class LentiviralRiskAssessor:
    """慢病毒包装风险评估器"""
    
    def __init__(self):
        self.risk_keywords = {
            "high": {
                "antiviral": ["interferon", "ifn", "antiviral", "innate immunity", "rig-i", "tlr", "sting", "mavs", "irf"],
                "toxic": ["lethal", "essential", "cell death", "apoptosis", "toxic", "fatal"],
                "proliferation": ["cell cycle arrest", "growth inhibition", "anti-proliferative", "tumor suppressor", "contact inhibition"],
                "structure": ["transmembrane domain", "secreted protein", "extracellular matrix", "collagen"]
            },
            "medium": {
                "signaling": ["kinase", "phosphatase", "signal transduction", "pathway"],
                "transcription": ["transcription factor", "nuclear receptor", "epigenetic"]
            }
        }
    
    def assess_by_function(self, gene_description: str, phenotype: str) -> Dict:
        """根据基因功能评估风险"""
        desc_lower = gene_description.lower()
        risks = []
        risk_level = "low"
        
        for level, categories in self.risk_keywords.items():
            for category, keywords in categories.items():
                matched = [kw for kw in keywords if kw in desc_lower]
                if matched:
                    risks.append(f"{category}: 检测到关键词 {matched[:2]}...")
                    if level == "high":
                        risk_level = "high"
                    elif level == "medium" and risk_level != "high":
                        risk_level = "medium"
        
        if "必需" in phenotype or "lethal" in phenotype.lower() or "essential" in phenotype.lower():
            risks.append("致死性: 必需基因，敲除可能导致细胞死亡")
            risk_level = "high"
        
        return {
            "risk_level": risk_level,
            "risks": risks,
            "recommendation": self._get_recommendation(risk_level, risks)
        }
    
    def _get_recommendation(self, risk_level: str, risks: List[str]) -> str:
        """根据风险等级给出建议"""
        if risk_level == "high":
            return "⚠️ 高风险：建议使用诱导型系统(Tet-On/Off)或暂时性敲低(shRNA)，避免直接KO"
        elif risk_level == "medium":
            return "⚡ 中等风险：可进行KO但需密切监测细胞状态，建议准备诱导型备选方案"
        else:
            return "✅ 低风险：标准KO方案适用，预期可获得稳定细胞系"
    
    def assess_by_literature(self, literature_data: Dict) -> Dict:
        """根据文献评估包装可行性"""
        oe_count = literature_data.get("overexpression", {}).get("count", 0)
        kd_count = literature_data.get("knockdown", {}).get("count", 0)
        ko_count = literature_data.get("knockout", {}).get("count", 0)
        
        evidence = {
            "overexpression": {"available": oe_count > 0, "count": oe_count, "method": "慢病毒包装"},
            "knockdown": {"available": kd_count > 0, "count": kd_count, "method": "shRNA/siRNA"},
            "knockout": {"available": ko_count > 0, "count": ko_count, "method": "CRISPR/Cas9"}
        }
        
        if oe_count > 10:
            packaging_feasibility = "high"
        elif oe_count > 0:
            packaging_feasibility = "medium"
        else:
            packaging_feasibility = "unknown"
        
        return {
            "evidence": evidence,
            "packaging_feasibility": packaging_feasibility,
            "has_precedent": oe_count > 0 or kd_count > 0 or ko_count > 0
        }
    
    def extract_sequences_from_literature(self, articles: List[Dict]) -> Dict:
        """从文献中提取shRNA/siRNA/sgRNA序列"""
        sequences = {
            "shrna": [],
            "sirna": [],
            "sgrna": []
        }
        
        patterns = {
            "shrna": r'(?:shRNA|shRNA\s+sequence)[:\s]+([ACGTU]{19,23})',
            "sirna": r'(?:siRNA|siRNA\s+sequence)[:\s]+([ACGTU]{19,23})',
            "sgrna": r'(?:sgRNA|gRNA|guide\s+RNA)[:\s]+([ACGTU]{20,23})',
            "target_seq": r'target\s+sequence[:\s]+([ACGTU]{19,23})',
            "forward": r'[Ff]orward[:\s]+([ACGTU]{19,23})',
            "sense": r'[Ss]ense[:\s]+([ACGTU]{19,23})'
        }
        
        for article in articles:
            text = article.get("title", "") + " " + article.get("abstract_snippet", "")
            
            for seq_type, pattern in patterns.items():
                matches = re.findall(pattern, text)
                for match in matches:
                    seq = match[0] if isinstance(match, tuple) else match
                    seq = seq.upper().replace("U", "T")
                    
                    if len(seq) >= 19 and len(seq) <= 23 and all(c in "ATCG" for c in seq):
                        entry = {
                            "sequence": seq,
                            "pmid": article.get("pmid"),
                            "title": article.get("title", "")[:50] + "..." if len(article.get("title", "")) > 50 else article.get("title", ""),
                            "type": seq_type
                        }
                        
                        if "shrna" in seq_type.lower():
                            if not any(e["sequence"] == seq for e in sequences["shrna"]):
                                sequences["shrna"].append(entry)
                        elif "sirna" in seq_type.lower():
                            if not any(e["sequence"] == seq for e in sequences["sirna"]):
                                sequences["sirna"].append(entry)
                        elif "grna" in seq_type.lower() or "sgrna" in seq_type.lower():
                            if not any(e["sequence"] == seq for e in sequences["sgrna"]):
                                sequences["sgrna"].append(entry)
        
        for seq_type in sequences:
            sequences[seq_type] = sequences[seq_type][:5]
        
        return sequences


# ==================== NCBI 数据获取模块 ====================

class BioDataFetcher:
    def __init__(self, email: str, api_key: str = ""):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.uniprot_base = "https://rest.uniprot.org/uniprotkb/search.json"
        self.headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'}
        self.addgene_scraper = AddgeneScraper()
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
                "description": summary[:800] if summary else "无描述",
                "phenotype": phenotype,
                "chromosome": gene_entry.get("Entrezgene_location", [{}])[0].get("Gene-location", {}).get("Gene-location_chromosome", "N/A"),
                "status": "success"
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}
    
    def get_uniprot_info(self, gene_symbol: str, species: str) -> Dict:
        """增强版 UniProt 查询，支持多种匹配策略"""
        try:
            species_map = {
                "Homo sapiens": ("human", 9606),
                "Mus musculus": ("mouse", 10090),
                "Rattus norvegicus": ("rat", 10116)
            }
            org_name, tax_id = species_map.get(species, (species.lower(), None))
            
            queries = [
                f"gene:{gene_symbol}+organism_id:{tax_id}",
                f"gene:{gene_symbol}+organism:{org_name}",
                f"{gene_symbol}+organism:{org_name}",
                gene_symbol
            ]
            
            for query in queries:
                try:
                    params = {
                        "query": query,
                        "fields": "accession,gene_names,length,cc_subcellular_location,sequence,protein_name",
                        "format": "json",
                        "size": 5
                    }
                    
                    response = requests.get(
                        self.uniprot_base, 
                        params=params, 
                        headers=self.headers, 
                        timeout=15
                    )
                    data = response.json()
                    
                    if data.get("results"):
                        best_match = None
                        gene_names = []
                        
                        for protein in data["results"]:
                            genes = protein.get("genes", [])
                            for g in genes:
                                if g.get("geneName"):
                                    gene_names.append(g["geneName"].get("value", "").upper())
                                if g.get("synonyms"):
                                    for syn in g["synonyms"]:
                                        gene_names.append(syn.get("value", "").upper())
                            
                            if gene_symbol.upper() in gene_names:
                                best_match = protein
                                break
                        
                        if not best_match:
                            best_match = data["results"][0]
                        
                        accession = best_match.get("primaryAccession", "")
                        seq_length = best_match.get("sequence", {}).get("length", 0)
                        
                        loc_text = ""
                        comments = best_match.get("comments", [])
                        for comment in comments:
                            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                                locations = comment.get("subcellularLocations", [])
                                locs = [loc.get("location", {}).get("value", "") 
                                       for loc in locations if loc.get("location")]
                                loc_text = "; ".join([l for l in locs if l])
                        
                        cds_length = seq_length * 3 if seq_length else 0
                        
                        return {
                            "uniprot_id": accession,
                            "protein_length": seq_length,
                            "cds_length_bp": cds_length,
                            "subcellular_location": loc_text or "未标注",
                            "protein_name": best_match.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
                            "status": "success",
                            "match_type": "exact" if gene_symbol.upper() in gene_names else "partial"
                        }
                except:
                    continue
            
            return {
                "status": "not_found",
                "error": f"UniProt 未找到 {gene_symbol}，请确认基因符号"
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
                raise ValueError("未配置 BAILIAN_API_KEY")
            
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
        self.risk_assessor = LentiviralRiskAssessor()
    
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
            
            hpa_gene_data = {}
            if species == "Homo sapiens":
                st.text("检索 HPA 蛋白表达数据...")
                hpa_gene_data = self.fetcher.hpa_data.get_gene_data(gene_symbol)
            time.sleep(0.5)
            
            st.text("检索文献...")
            literature = self._search_all_constructs(gene_symbol, cell_line)
            time.sleep(0.5)
            
            st.text("评估慢病毒风险...")
            lentiviral = self._assess_lentiviral_comprehensive(
                gene_symbol, ncbi_info, uniprot_info, literature
            )
            
            ai_analysis = {}
            try:
                ai_analyzer = AIAnalyzer()
                
                with st.spinner("🤖 通义千问正在深度分析..."):
                    ai_analysis["gene_assessment"] = ai_analyzer.assess_gene_essentiality(
                        ncbi_info, uniprot_info
                    )
                    
                    if cell_line and literature.get("specific_cell", {}).get("found"):
                        for ctype in ["overexpression", "knockdown", "knockout"]:
                            if literature["specific_cell"][ctype]["articles"]:
                                ai_result = ai_analyzer.analyze_literature_deep(
                                    literature["specific_cell"][ctype]["articles"], gene_symbol, ctype
                                )
                                literature["specific_cell"][ctype]["ai_analysis"] = ai_result
                    
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
    
    def _assess_lentiviral_comprehensive(self, gene_symbol: str, ncbi_info: Dict, 
                                        uniprot_info: Dict, literature: Dict) -> Dict:
        """综合评估慢病毒适用性"""
        cds_len = uniprot_info.get("cds_length_bp", 0)
        
        # 1. CDS长度评估
        if cds_len > 2500:
            cds_risk = {"level": "high", "suitable": False, "reason": f"CDS过长({cds_len}bp)，远超2000bp推荐值"}
        elif cds_len > 2000:
            cds_risk = {"level": "medium", "suitable": True, "reason": f"CDS较长({cds_len}bp)，建议使用第三代系统"}
        elif cds_len == 0:
            cds_risk = {"level": "unknown", "suitable": True, "reason": "无法获取CDS长度信息"}
        else:
            cds_risk = {"level": "low", "suitable": True, "reason": f"CDS长度合适({cds_len}bp)"}
        
        # 2. 基因功能风险评估
        function_risk = self.risk_assessor.assess_by_function(
            ncbi_info.get("description", ""),
            ncbi_info.get("phenotype", "")
        )
        
        # 3. 文献证据评估
        lit_evidence = self.risk_assessor.assess_by_literature(literature)
        
        # 4. 提取序列
        sequences = {}
        for ctype in ["knockdown", "knockout"]:
            if literature[ctype].get("articles"):
                seqs = self.risk_assessor.extract_sequences_from_literature(
                    literature[ctype]["articles"]
                )
                if any(seqs.values()):
                    sequences[ctype] = seqs
        
        # 综合建议
        recommendations = []
        warnings = []
        
        if cds_risk["level"] == "high":
            warnings.append(f"⚠️ 长度风险: {cds_risk['reason']}")
            recommendations.append("建议使用split-vector系统或选择其他递送方式（如转座子）")
        elif cds_risk["level"] == "medium":
            warnings.append(f"⚡ 长度警告: {cds_risk['reason']}")
        
        if function_risk["risk_level"] == "high":
            warnings.append(f"🚨 功能风险: {', '.join(function_risk['risks'][:2])}")
            recommendations.append("强烈建议使用诱导型表达系统（Tet-On/Off）")
        elif function_risk["risk_level"] == "medium":
            warnings.append(f"⚡ 功能注意: {', '.join(function_risk['risks'][:1])}")
        
        if not lit_evidence["has_precedent"]:
            warnings.append("📚 文献缺乏: 未找到该基因的病毒包装文献记录")
            recommendations.append("建议先进行小规模包装测试")
        else:
            if lit_evidence["evidence"]["overexpression"]["available"]:
                recommendations.append("✅ 文献支持: 已有成功过表达记录")
        
        # 总体评级
        if cds_risk["level"] == "high" or function_risk["risk_level"] == "high":
            overall_rating = "❌ 高风险（不推荐标准方案）"
        elif cds_risk["level"] == "medium" or function_risk["risk_level"] == "medium":
            overall_rating = "⚠️ 中等风险（需优化方案）"
        else:
            overall_rating = "✅ 低风险（标准方案适用）"
        
        return {
            "cds_assessment": cds_risk,
            "function_risk": function_risk,
            "literature_evidence": lit_evidence,
            "sequences": sequences,
            "warnings": warnings,
            "recommendations": recommendations,
            "overall_rating": overall_rating,
            "overall_suitable": cds_risk["suitable"] and function_risk["risk_level"] != "high"
        }
    
    def _search_all_constructs(self, gene_symbol: str, cell_line: Optional[str]) -> Dict:
        """双分支文献检索"""
        results = {}
        
        # 分支1：特定细胞系
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
                            
                            methods = [kw for kw in ["lentiviral", "transfection", "electroporation", "transduction"] if kw in full_text]
                            
                            article_info = {
                                "pmid": pmid,
                                "title": title_display,
                                "methods": ", ".join(methods) if methods else "未明确提及",
                                "cell_line": cell_line,
                                "abstract_snippet": abstract[:300] if len(abstract) > 300 else abstract
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
                        "message": f"在 {cell_line} 中找到 {len(pmids)} 篇文献"
                    }
                else:
                    results["specific_cell"] = {"found": False, "message": f"未在 {cell_line} 中找到相关研究"}
            except Exception as e:
                results["specific_cell"] = {"found": False, "message": f"检索失败: {e}"}
        else:
            results["specific_cell"] = {"found": False, "message": "未输入细胞系名称"}
        
        # 分支2：通用文献
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
            "lentiviral_suitable": lentiviral.get("overall_suitable"),
            "lentiviral_risk": lentiviral.get("function_risk", {}).get("risk_level"),
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
    st.markdown("**2000bp包装极限版 + 功能风险评估 + 序列提取**")
    
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
        
        hpa_checker = HPAGeneData()
        if hpa_checker.check_data_available():
            st.success("✅ HPA 人类蛋白数据已加载")
        else:
            st.warning("⚠️ HPA 数据未配置")
            st.caption("首次使用将自动下载（约30MB）")
        
        st.info("""
        **系统功能：**
        - ✅ NCBI/UniProt/HPA 多源数据
        - ✅ CDS长度 + 功能风险双评估
        - ✅ 文献包装证据检索
        - ✅ shRNA/siRNA/sgRNA序列提取
        - ✅ Addgene质粒查询
        - ✅ AI智能分析
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
            st.exception(e)
    
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
    
    tabs = st.tabs(["🧬 基因功能", "🦠 慢病毒风险评估", "📚 文献与序列", "🧪 实验资源"])
    
    # Tab 1: 基因功能
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
                if prot_data.get("protein_name"):
                    st.write(f"**蛋白名称:** {prot_data.get('protein_name')[:50]}...")
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
                if hpa_data.get('rna_tissue_specificity'):
                    rna = str(hpa_data['rna_tissue_specificity'])
                    if len(rna) > 50:
                        rna = rna[:47] + "..."
                    st.write(f"**RNA表达:** {rna}")
                with st.expander("亚细胞定位"):
                    st.write(hpa_data.get('subcellular_location', '未标注'))
                st.caption(f"[HPA详情页]({hpa_data.get('hpa_link', '#')})")
            else:
                if info['species'] == "Homo sapiens":
                    st.info("未找到HPA数据")
                else:
                    st.info("HPA仅支持人类基因")
    
    # Tab 2: 慢病毒风险评估
    with tabs[1]:
        lv = result["lentiviral_assessment"]
        
        st.subheader("综合风险评估")
        
        col1, col2, col3 = st.columns(3)
        col1.metric("CDS长度风险", lv['cds_assessment']['level'].upper(), 
                   help=lv['cds_assessment']['reason'])
        col2.metric("功能风险等级", lv['function_risk']['risk_level'].upper(),
                   help="基于基因功能描述的风险评估")
        col3.metric("总体评级", lv['overall_rating'].split()[0])
        
        # 警告信息
        if lv['warnings']:
            st.error("**⚠️ 风险提示:**")
            for warning in lv['warnings']:
                st.write(f"- {warning}")
        
        # 建议
        if lv['recommendations']:
            st.success("**💡 专家建议:**")
            for rec in lv['recommendations']:
                st.write(f"- {rec}")
        
        # 功能风险详情
        with st.expander("查看功能风险详情"):
            if lv['function_risk']['risks']:
                for risk in lv['function_risk']['risks']:
                    st.write(f"- {risk}")
            else:
                st.write("未检测到特殊功能风险")
            st.info(f"**策略建议:** {lv['function_risk']['recommendation']}")
        
        # 文献证据
        with st.expander("查看文献包装证据"):
            ev = lv['literature_evidence']['evidence']
            for construct, data in ev.items():
                status = "✅" if data['available'] else "❌"
                st.write(f"{status} **{construct}**: {data['count']}篇文献 ({data['method']})")
    
    # Tab 3: 文献与序列
    with tabs[2]:
        literature = result["cell_line_constructs"]
        lv = result["lentiviral_assessment"]
        
        # 文献列表
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
        
        # 提取的序列
        if lv.get('sequences'):
            st.divider()
            st.subheader("🧬 文献报道的靶点序列")
            
            for construct_type, seqs in lv['sequences'].items():
                if any(seqs.values()):
                    with st.expander(f"{construct_type.upper()} 序列", expanded=True):
                        for seq_type, entries in seqs.items():
                            if entries:
                                st.text(f"{seq_type.upper()} 序列 (来自文献):")
                                for entry in entries:
                                    cols = st.columns([2, 1, 3])
                                    cols[0].code(entry['sequence'])
                                    cols[1].caption(f"PMID: {entry['pmid']}")
                                    cols[2].caption(entry['title'][:40])
        
        # 通用文献统计
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
